using TaylorModels
include("flowstar.jl")

"""Convert the given box to a TM initial set.

    A box also represents a TM initial set. This is a different,
    but equivalent representation to a TM initial set.
"""
function tm_initial_set(box::IntervalBox)
    # Each variable should have a matching box element.
    @assert(get_numvars() == length(box))

    # Do the conversion: x  =>  (x, [0, 0])  given the domain of x.
    zd = zero(box)
    vars = get_variables()
    # FIXME: This also adds an entry for variable t, but it should not?
    #       Because an ODE with components x,y should have a 2-dimensional
    #       init set Xi, with components only for x,y and not t?
    return map((var) -> TaylorModelN(var, 0..0, zd, box), vars)
end


"""Given a single TM (p, I), compute the corresponding flowpipe.

    Given an initial set Xi, the flowpipe is computed as (p(Xi, t), I).
"""
function flowpipe(p::TaylorN, I::Interval, Xi::Vector, domain::IntervalBox)
    Xi = copy(Xi)
    # Compute the flowpipe Fi = (p(Xi), I) = (ps, Is + I)
    #             where p(Xi) = (ps, Is)
    tm = p(Xi)
    ps = polynomial(tm)
    Is = remainder(tm)
    zd = zero(domain)
    return TaylorModelN(ps, Is + I, zd, domain)
end


"""Compute the next iteration's initial set Xi based on the
    current flowpipe Fi.
"""
function initial_set(Fi::TaylorModelN, δi::Float64, domain::IntervalBox)
    # FIXME: AFAIK the 'evaluate' function for TaylorN of TaylorSeries
    #       does not allow mixing value types. For example:
    # x, y, a, b = set_variables("x y a b", order=3)
    # p = TaylorN(x + y)
    # p([1, 2, 3, 4])   # This works fine.
    # p([a, b, x, y])   # This works fine too.
    # p([1, b, 3, y])   # ERROR.
    # FIXME: The solution? A dirty hack, hooray.

    # Construct a constant TaylorN expression of δi.
    δi_cte = TaylorN(δi, get_order())
    # To compute Xi, fix t=δi in Fi.
    # FIXME: Assume t is the last variable and exclude it from the list.
    vars_except_t = get_variables()[1:end-1]
    values = vcat(vars_except_t, [δi_cte])
    println("subst values = $values ($(typeof(values)))")
    ps = polynomial(Fi)(values)
    Is = remainder(Fi)
    zd = zero(domain)
    return TaylorModelN(ps, Is, zd, domain)
end

"""The standard TM integration algorithm.

    Compute an overapproximation of the true flow of the
    system of ODEs represented by the given vector field,
    as a sequence of flowpipes over partial time horizons
    [0, δi] of the full time horizon [0, Δ].
"""
function tm_integration(f, domain, k::Integer, J, Δ::Float64,
                        TIME_STEP_SIZE::Float64,
                        NR_CONTRACTIVENESS_TRIES::Integer,
                        NR_REFINEMENTS::Integer,
                        SCALE::Float64)
    @assert(Δ > 0)  # The time horizon must not be [0, 0].
    @assert(TIME_STEP_SIZE > 0)
    @assert(NR_CONTRACTIVENESS_TRIES >= 0)
    @assert(NR_REFINEMENTS >= 0)

    # Derive the corresponding TM initial set from the interval initial set.
    X0 = tm_initial_set(domain)
    zd = zero(domain)

    # Start the TM integration loop.
    remaining_time::Float64 = Δ
    initial_sets::Array{Any} = [X0]
    Fi = nothing    # The current flowpipe.
    Xi = X0         # The current initial set.
    # FIXME: Would a for-loop be cleaner?
    while remaining_time > 0.0
        δi = min(remaining_time, TIME_STEP_SIZE)
        remaining_time -= δi
        println("time step [0, $δi] (remaining Δ: $remaining_time)")

        # Step 1: generate the poly approximation of the flow.
        p = tay_poly(f, k)

        # Step2: generate a safe, refined remainder interval.
        I = tay_model_error(f, p, domain, k, J,
                            NR_CONTRACTIVENESS_TRIES,
                            NR_REFINEMENTS,
                            SCALE)

        # get a copy of the t variable before messing up order
        t = get_variables()[end]

        # Step 3: Compute the flowpipe and next initial set.
        # To account for composition, we'll up the current order
        old_order = get_order()
        hgr_order = old_order * maximum(get_order.(p))
        vars = set_variables(get_variable_string(), order=hgr_order)

        println("p = $p")
        println("I = $I")

        # TODO: Compute the flowpipe components.
        # FIXME: Does TaylorModels actually add the enclosures of truncated
        # terms during TM arithmetic to the remainders? Because the
        # substitution below requires this!
        tmv = zip(p, I)
        Fi = map(((pj, Ij),) -> flowpipe(pj, Ij, Xi, domain), tmv)
        println("Fi = $Fi")

        # TODO: Compute the initial set components.
        Xi = map((fpi) -> initial_set(fpi, δi, domain), Fi)
        
        dummy_t_tm = TaylorModelN(t, 0..0, zd, domain)
        append!(Xi, dummy_t_tm)   # FIXME: Append dummy

        println("Xi = $Xi")
        append!(initial_sets, Xi)

        # We now restore the order.
        vars = set_variables(get_variable_string(), order=old_order)
        # break   # FIXME: REMOVE, this is here for testing purposes (i.e. faster testing by breaking loop)
    end

    return initial_sets
end


if abspath(PROGRAM_FILE) == @__FILE__

    # Full TM integration, based on Example 3.3.2 and Example 3.3.6
    k = 3   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 5
    NR_REFINEMENTS           = 1
    SCALE                    = 2.0
    TIME_STEP_SIZE = 0.02
    Δ = 0.2     # The finite time horizon
    vars = set_variables("x y t", order=k)
    # The vector field f of the ODEs:
    #   f[1] = 1 + y
    #   f[2] = -x^2
    f = [1 + vars[2],  # x
        -vars[1]^2]   # y
    domain = IntervalBox([-1..1,      # x
                -0.5..0.5,  # y
                0..0.02])   # t
    # Initial remainder estimate J, a hyperrectangle
    J = fill(-0.1..0.1, length(f))

    initial_sets = tm_integration(f, domain, k, J, Δ,
                               TIME_STEP_SIZE,
                               NR_CONTRACTIVENESS_TRIES,
                               NR_REFINEMENTS,
                               SCALE)
    println("init sets = $initial_sets")


    # x, y, a, b = set_variables("x y a b", order=3)
    # p = TaylorN(x + y)
    # println("p([1, 2, 3, 4]): $(p([1, 2, 3, 4]))")
    # println("p([a, b, x, y]): $(p([a, b, x, y]))")
    # println("p([1, b, 3, y]): $(p([1, b, 3, y]))")


    # X0 = tm_initial_set(domain)
    # println(X0)



    # p = tay_poly(f, k)

    #     # To account for composition, we'll up the current order
    #     old_order = get_order()
    #     hgr_order = old_order * maximum(get_order.(p))
    #     vars = set_variables(get_variable_string(), order=hgr_order)

    # println(flowpipe(p[1], -1..1, X0, domain))

    #     # We now restore the order.
    #     vars = set_variables(get_variable_string(), order=old_order)

end
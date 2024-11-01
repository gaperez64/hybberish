using TaylorSeries  # IntervalArithmetic
using TaylorModels

# Assumes the last variable is t
"""Generate the Taylor polynomial approximation part of the i-th flowpipe
    via Lie derivatives.
"""
function tay_poly(f, k)
    vars = get_variables()
    t = vars[end]
    # Let's pad f with a 1 at the end for t
    fp1 = copy(f)
    push!(fp1, 1)
    # Also, prepare a valuation vector with t=0
    # FIXME: We're cheating to obtain a zero with same order as the other
    # variables
    val0 = copy(vars)
    val0[end] = t - t
    # Prepare a first lie derivative and the result
    g = copy(vars)
    deleteat!(g, length(vars))
    # Start a vector function for the result
    res = copy(g)
    
    for i = 1:k
        g = TaylorSeries.jacobian(g, vars) * fp1
        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
        res += term
    end
    
    return res
end


"""Get the TaylorSeries variable names as a space-separated string."""
function get_variable_string()
    return join(get_variable_names(), " ")
end


"""Compute the remainder of the TM extension of the picard operator:
    Pf((p, I)) = J.
"""
function picard_tm_extension(f, tmv, domain, k)
    # FIXME: Is a copy needed? The substitution f(g(x, t), t)
    # seems to have side effects? Or is it just getting late :/
    tmv = copy(tmv)

    # Perform the substitution operation of the Picard operator:
    #   f(g(x, s), s).
    ftm = map((fj) -> fj(tmv), f)

    # Compute the integral's remainder part:
    #   (Int(pe) + I) * [0, t].
    # FIXME: Is it fine to assume the t interval is the last one in the box?
    tdom = domain.v[end]
    # FIXME: Is readable at all?
    Int(p_error_terms, dom) = sum(perr(dom) for perr in p_error_terms)
    pe(tm::TaylorModelN) = polynomial(tm)[k:end]
    integrated = map((tmj) -> (Int(pe(tmj), domain) + remainder(tmj)) * tdom, ftm)

    return integrated
end


"""A broadcasting of the TM remainder function across a TM vector."""
remainders(tmv) = remainder.(tmv)


"""Compute a safe remainder interval for the i-th flowpipe.

    The TM extension of the Picard operator, Pf, can be used to
    refine the remainder estimate of the i-th flowpipe.
    Suppose we compute the TM (p, J_{0}) before refinement.
    This function computes a safe remainder by first finding a
    contractive remainder
        (p, J_{m})  where m >= 0
    and then continually applying refinement
        Pf((pl, J_{j})) = (pl, J_{j+1})  for n = m, ..., n
    to finally obtain the contractive, refined TM
        (p, J_{n}).
    Note that m, n are usually not known beforehand, but instead
    become known when the contractiveness and refinement iteration
    converge or fail.
"""
function tay_model_error(f, p, domain, k::Integer, J,
                         NR_CONTRACTIVENESS_TRIES::Integer,
                         MAX_WIDTH_THRESHOLD::Float64,
                         NR_REFINEMENTS::Integer,
                         SCALE::Float64)
    zd = zero(domain)

    # To account for composition, we'll square the current order
    old_order = get_order()
    # TODO: This is not exactly squaring since you multiply
    # two potentially different values?
    sqd_order = old_order * maximum(get_order.(p))
    vars = set_variables(get_variable_string(), order=sqd_order)
    # FIXME: The time variable MUST be of the old variable order.
    # We define a dummy TM with it, and else TaylorModels'
    # TM multiplication throws an assertion error about
    # the terms having too large an order.
    t = get_variables(old_order)[end]
    # t = vars[end] # FIXME: Uncomment to recreate the error described above.

    # TODO: start loop
    #   ==> [Thomas]; Done?

    # Construct the TM vector now.
    # Note: Anonymous function tuple destructuring has unique syntax:
    #   https://discourse.julialang.org/t/argument-destructuring-and-anonymous-functions/24893
    compose_tm = ((poly, rem),) -> TaylorModelN(poly, rem, zd, domain)
    tmv = map(compose_tm, zip(p, J))

    # Add a dummy TM (t, [0, 0]) for the time variable because the substitution
    # requires one TM for every variable, which includes the time variable.
    # FIXME: Can this ugliness be circumvented?
    dummy_t_tm = TaylorModelN(t, 0..0, zd, domain)
    append!(tmv, dummy_t_tm)

    # TODO: truncate and compute error
    #   ==> [Thomas]; NOT done I think!
    # TODO: end loop after contraction check
    #   ==> [Thomas]; Done? See "break" in for-loop below.
    J0 = nothing
    J1 = nothing
    for _ in 1:NR_CONTRACTIVENESS_TRIES
        J0 = remainders(tmv)[1:end-1]   # FIXME: drop the dummy interval of t

        # Width([I1, ..., In]) = max(Width(I1), ..., Width(In)).
        # If any remainder exceeds the width threshold, then return FAIL.
        if any(diam(j0) >= MAX_WIDTH_THRESHOLD for j0 in J0)
            @assert(false)
        end

        # Perform the first refinement, which to next test for contractiveness.
        J1 = picard_tm_extension(f, tmv, domain, k)

        # Test contractiveness.
        if all(issubset.(J1, J0))
            # Contractiveness test succeeded, pass along the contractive (safe)
            # remainder J1.
            tmv = map(compose_tm, zip(p, J1))
            append!(tmv, dummy_t_tm)    # FIXME: Still appending dummy TM :(
            break
        end

        # Contractiveness test failed, widen all initial remainders.
        tmv = map(compose_tm, zip(p, remainders(tmv)*SCALE))
        append!(tmv, dummy_t_tm)    # FIXME: Still appending dummy TM :(
    end

    # TODO: add loop to continue contracting for a while
    #   ==> [Thomas]; Done?
    for _ in 1:NR_REFINEMENTS
        Jn = picard_tm_extension(f, tmv, domain, k)
        tmv = map(compose_tm, zip(p, Jn))
        append!(tmv, dummy_t_tm)    # FIXME: Still appending dummy TM :(
    end

    # Collect the results.
    tmv = tmv[1:end-1]    # FIXME: Drop dummy TM for final output.
    errors = remainders(tmv)

    # We now restore the order.
    vars = set_variables(get_variable_string(), order=old_order)
    return errors
end

# TODO: For now, we require a custom version of the TaylorModels package
#   https://github.com/gaperez64/TaylorModels.jl

#TODO: Clean below and call the above
#
# Example 3.3.6
k = 3   # The TM arithmetic and truncation order
NR_CONTRACTIVENESS_TRIES = 5
MAX_WIDTH_THRESHOLD      = 10.0
NR_REFINEMENTS           = 1
SCALE                    = 2.0
vars = set_variables("x y t", order=k)
# The vector field f of the ODEs:
#   f[1] = 1 + y
#   f[2] = -x^2
f = [1 + vars[2], -vars[1]^2]
# The polynomial approx:
#   p[1] = x + t + yt
#   p[2] = y - (x^2)t - xt^2 - (1/3)t^3
p = [
    vars[1] + vars[3] + vars[2] * vars[3],
    vars[2] - vars[1]^2 * vars[3] - vars[1] * vars[3]^2 - 1/3 * vars[3]^3
]
domain = IntervalBox([-1..1,-0.5..0.5,0..0.02])
# TODO: [Thomas]: I personally prefer the initial remainder estimate J (sometimes
# also called J0) to be an input to the TM integration function. It can be fixed
# as a constant in some config file, but preferably not in the code itself?
J = fill(-0.1..0.1, length(p))
I = tay_model_error(f, p, domain, k, J,
                    NR_CONTRACTIVENESS_TRIES,
                    MAX_WIDTH_THRESHOLD,
                    NR_REFINEMENTS,
                    SCALE)
println(" safe remainders = $I")


# Example 3.3.2
vars = set_variables("x y t", order=4)
f = [1 + vars[2], -vars[1]^2]
p = tay_poly(f, 4)
println("poly from new TM = $p")

include("taylor_models/BasicTaylorModels.jl")


"""Generate the Taylor polynomial approximation of the true flow specified by
   given dynamics (ODEs) up to the given degree via Lie derivatives.

   We require TaylorN representations of the variables to be supplied by the
   caller. This allows the caller to decide their order; the variables' orders
   influence the orders of the output polynomials.

   @param[in]    f The Taylorized dynamics (ODEs).
   @param[in]    k The degree (order) of the resulting polynomials.
   @param[in] vars The TaylorN objects representing the variables.
                   Assume the last element represents the time variable.
   @return The polynomial approximations of the true flow.
"""
function tay_poly(f::Vector{TaylorN{N}}, k::Integer, vars::Vector{TaylorN{N}}) where {N <: Number}
    vars_orders = get_order.(vars)
    # The function `evaluate(::TaylorN, ::Vector{TaylorN})` does not play nice
    # with different orders for the substitution values. So, require all orders
    # to be the same.
    @assert all( vars_orders[1] .== vars_orders )
    # Taylor series arithmetic propagates the lowest order of its operands.
    # To generate order k polynomials, all variables must be at least order k.
    @assert all( k .<= vars_orders )
    # The variables contain an additional last component: the time variable.
    @assert length(f) == length(vars)-1

    t = vars[end]
    # Let's pad f with a 1 at the end for t
    fp1 = copy(f)
    push!(fp1, 1)
    # Also, prepare a valuation vector with t=0
    # NOTE: We're cheating to obtain a zero with same order as the other
    # variables
    val0 = copy(vars)
    val0[end] = t - t
    # Prepare a first lie derivative and the result,
    # in particular we remove the time variable t
    g = copy(vars)
    deleteat!(g, length(vars))
    # Start a vector function for the result
    res = copy(g)
    for i = 1:k
        g = TaylorSeries.jacobian(g, vars) * fp1
        println("Lie derivative:")
        println(map((h) -> evaluate(h, vars), g))
        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
        res += term
    end
    return res
end


"""Compute the remainder of the TM extension of the picard operator.

    The term "TM extension" refers to the fact that all mathematical operators
    in the Picard operator are replaced with their TM arithmetic counterparts.
    This implies that the vector field f and function g w.r.t. which the normal
    Picard operator is computed, must be TMs instead of functions when working
    with the TM extension of the Picard operator P.
    Suppose

        P_f(g) = x * int_0^t f(g(x, t), t) dt

    is the normal Picard operator. Then

        P_F(G)
        = x * int_0^t F(G) dt
        = x * int_0^t (p, I) \\circ (q, J) dt
        = x * int_0^t p(q, J) + I dt
        = x * int_0^t (r, K) + I dt
        = x * int_0^t (r, K + I) dt
        = (p, K + I)

    is the TM extension of the Picard operator, where F = (p, I) is a TM
    (Taylorization plus Lagrange remainder) of the vector field f and
    G = (q, J) is the TM representing the function g.

    @param[in] vector_field_tms The TM representation F of the vector field f.
    @param[in]     function_tms The TM representation G of the function g.
    @param[in]           domain The domains of all ODE variables, including time.
    @param[in]                k The Taylor polynomial and truncation order.
    @return The remainders of the TM result of the TM extension of the Picard
     operator; the interval vector K.
"""
function picard_tm_extension(vector_field_tms::Vector{T}, function_tms, domain, k) where T <: TaylorModelN
    #= Step (1), setup. =#
    ode_polynomials::Vector{TaylorN} = polynomial.(vector_field_tms)
    ode_remainders::Vector{Interval} = remainder.(vector_field_tms)
    # FIXME: Is a copy needed? f(g(x,t),t) has side effects?
    function_tms = copy(function_tms)

    #= Step (2), perform the substitution operation part of the
      TM extension of Picard operator:
        (p, I) \circ (q, J) = p(q, J) + I = (r, K) + I = (r, K + I)
      The normal Picard operator expresses this as f(g(x, t), t).
    =#
    #= Step (2.1), perform p(q, J) = (r, K) =#
    substitution_tms = map((fj) -> evaluate(fj, function_tms), ode_polynomials)
    #= Step (2.2), compute (r, K) + I = r + K + I = (r, K + I)
      Add the error of aproximating the true vector field. Here we have NOT
      yet reached where we can use formula Xin Chen (p42, antiderivative formula)
      we are still just accounting for errors coming from the dynamics having
      been approximated by polynomials. Here we're instead using the formula
      for the error of TM composition:
        (p, I) \circ (q, J)
        = p(q, J) + I
        = (r, K) + I
        = (r, K + I)
    =#
    substitution_tms = map(
        ((tm, err),) -> TaylorModelN(
            polynomial(tm),      # r
            remainder(tm) + err, # K + I
            domain),
        zip(substitution_tms, ode_remainders))

    #= Step (3), apply the antiderivative formula. =#
    # FIXME: Once again we assume the t interval is the last one
    tdom = domain.v[end]
    # Prepare two helper functions to compute the antiderivative remainder.
    int_enclosure(poly_error_terms, domains) =
        sum(evaluate(perr, domains) for perr in poly_error_terms)
    pe(tm::TaylorModelN) = polynomial(tm)[k:end]
    #= Compute the integral's remainder part:
        (int_enclosure(pe) + I) * [t, t+\delta].
      where pe represents the sum of all truncated (error)
      terms of the polynomial integral result.
    =#
    return map(
        (tmj) -> (int_enclosure(pe(tmj), domain) + remainder(tmj)) * tdom,
        substitution_tms)
end


# Dynamics
f_dot(y, t) = -y - sin(t) + cos(t)

# Integration task specification
# a. Take delta_t = 1
# b. Construct a flowpipe consisting of 4 Taylor models, including the initial
#    one
# c. Work with order 4 (truncation degree / polynomial degree).
ord = 4
# Double the truncation degree to obtain the TaylorSeries max order.
# This accounts for order-related assertions applicable to Taylor
# series arithmetic.
ord_max = 2*ord

# Taylor variables (from TaylorSeries library)
set_variables("y t", order=ord_max)
vars = get_variables(ord)
y, t = vars

# Initial state variable bounds and domain
# y(0) = [1, 1]
# t(0) = [0, 0]
vals = IntervalBox(interval(1), interval(0))  # This is D_i in the maths

tstep = 0.1  # The fixed time step size.
scale = 2  # The scale factor for when contractiveness fails.


for _ = 0:10
    # Step 0: Taylorize the dynamics
    # We want to have a polynomial approximation of the dynamics centered around
    # the midpoint of the current values.
    fpoly = f_dot(y, t)
    println("taylorized vector field/dynamics:")
    println(fpoly)
    
    # Step 1: Obtain the polynomial part of the Taylor model
    p::Vector{TaylorN} = tay_poly([fpoly], ord, vars)
    println("polynomial part of TM:")
    println(p)
    
    # Step 2: Obtain the remainder/error interval of the TM
    rems = nothing
    tdom = vals[2].lo..(vals[2].lo+tstep) # [t0, t0+δ]
    doms = IntervalBox(vals[1], tdom)
    while true
        # Start Picard iteration, we need the candidate/guessed TM
        remainder_estimate = -0.1..0.1
        # based on the estimate, we want the flowpipe to be used as the domain for
        # the taylorization of the dynamics
        # FIXME: assuming the domain of time is a degenerate interval
        # NOTE: `fpipe` is F_i in the maths.
        fpipe = IntervalBox([
            p_i(doms) + remainder_estimate
            for p_i in p
        ]..., tdom)
        ytm = TaylorModelN(y, # polynomial
                interval(0),  # error remainder
                fpipe # domain hypercube
               )
        ttm = TaylorModelN(t, # polynomial
                interval(0),  # error remainder
                fpipe # domain hypercube
               )
        ftm = f_dot(ytm, ttm)
        println("poly version of dynamics, now with error")
        println(ftm)

        candidate_ttm = TaylorModelN(t, interval(0), doms)
        candidate_tm = TaylorModelN(p[1], remainder_estimate, doms)
        # Then we take the TM extension of the approx'd vector field composed
        # with the candidate TM.
        rems = picard_tm_extension([ftm], [candidate_tm, candidate_ttm], doms, ord)


        break # FIXME: Delete this! This allows non-contractive remainders!


        println("error interval part of TM:")
        if all(issubset.(rems, [remainder(candidate_tm)]))
            print(rems)
            print(" SUBSET ")
            println(remainder_estimate)
            println("Contractive!")
            break
        else
            print(rems)
            print(" NOT SUBSET ")
            println(remainder_estimate)
        end

        remainder_estimate *= 2
    end

    # FIXME: plot fpipe at this point
    #
    # Step 3: Get the new local values (and interval box) and update domain for next step
    # i.e. just change the domain of the time variable in doms
    valid_tm = TaylorModelN(p[1], rems[1], doms)
    doms = IntervalBox(doms[1], doms[2].hi..doms[2].hi)
    println("Full valid tm:")
    println(valid_tm)

    # NOTE: We are NOT evaluating valid_tm on vals because it complains about
    # it not being in the centered domain. Instead we manually compute the new
    # vals based on the polynomial part of valid_tm and its remainder.
    global vals = IntervalBox(valid_tm(doms)..., doms[2])
    println("                     vals: ", vals)
    println("\n=================================\n")
end

using TaylorSeries  # Uses IntervalArithmetic
using TaylorModels  # Version req'd: https://github.com/gaperez64/taylormodels.jl


# NOTE: Assumes the last variable is t
"""Generate the Taylor polynomial approximation part of a Taylor model for the
   given function and up to the given degree via Lie derivatives.
"""
function tay_poly(f::Vector{TaylorN{N}}, k::Integer) where {N <: Number}
    vars = get_variables()
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
    # FIXME: Hack to allow for higher degree terms in the intermediate computation
    y, t = set_variables("y t", order=k*2)
    substitution_tms = map((fj) -> evaluate(fj, function_tms), ode_polynomials)
    y, t = set_variables("y t", order=k)
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
            tm.x0, # Retain the expansion point.
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
        (int_enclosure(pe) + I) * [0, t].
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
# c. Work with order 4
ord = 4
domy = -3..3 # FIXME: Widened this from -2..2 because iscontained() failed in step 3.

# Taylor variables (from TaylorSeries library)
y, t = set_variables("y t", order=ord)

# Initial state variable bounds and domain
# y(0) = [1, 1]
# t(0) = [0, 0]
vals = IntervalBox(interval(1), interval(0))
doms = IntervalBox(domy, 0..0.5)

# FIXME: The very first time step size is governed by the time interval component
#  of the doms, NOT by the 'tstep' variable.
tstep = doms[2].hi - doms[2].lo # The fixed time step size.
scale = 2   # The scale factor for when contractiveness fails.

# FIXME: The second integration NEVER finds a contractive remainder for tstep=0.5!
#  So, manually make it smaller for testing purposes.
# tstep = 0.1

for it = 0:2
	# Step 0: Taylorize the dynamics
	# We want to have a polynomial approximation of the dynamics centered around
	# the midpoint of the current values.
	ytm = TaylorModelN(y, interval(0), vals, doms)
	println("ytm = ")
	println(ytm)
	ttm = TaylorModelN(t, interval(0), vals, doms)
	# FIXME: Hack to allow for higher degree terms in the intermediate computation
	set_variables("y t", order=ord*2)
	ftm = f_dot(ytm, ttm)
	# FIXME: We go back to the lower degree afterwards
	set_variables("y t", order=ord)
	println("taylorized vector field/dynamics:")
	println(ftm)

	# Step 1: Obtain the polynomial part of the Taylor model
	p = tay_poly([polynomial(ftm)], ord)
	println("polynomial part of TM:")
	println(p)

	# Step 2: Obtain the remainder/error interval of the TM
	# Start Picard iteration, we need the candidate/guessed TM
    remainder_estimate = -0.1..0.1
    rems = nothing
    while true
        candidate_tm = TaylorModelN(p[1], remainder_estimate, vals, doms)
        # Then we take the TM extension of the approx'd vector field composed
        # with the candidate TM
        rems = picard_tm_extension([ftm], [candidate_tm, ttm], doms, ord)
        break
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

	# Step 3: Get the new local values (and interval box) and update domain for next step
	# i.e. just change the domain of the time variable in doms
	valid_tm = TaylorModelN(p[1], rems[1], vals, doms)
    tdom = doms[2] # The time domain

    #= FIXME: The `evaluate(TaylorModelN, IntervalBox)` call asserts
      `iscontained(Interval, TaylorModelN)` which is defined as

          iscontained(a::Interval, tm::$TM) = a ⊆ centered_dom(tm)

      Then `centered_dom(TaylorModelN)` is defined as

          @inline centered_dom(tm::TaylorModelN) = domain(tm) .- expansion_point(tm)
    =#
    # FIXME: vals ⊆ centered_dom(valid_tm) MUST hold!
    println("                     vals: ", vals)
    println("         domain(valid_tm): ", domain(valid_tm))
    println("expansion_point(valid_tm): ", expansion_point(valid_tm))
    println("   centered_dom(valid_tm): ", centered_dom(valid_tm))

	global vals = IntervalBox(evaluate(valid_tm, vals), interval(tdom.hi))
    # FIXME: The y-component of vals keeps moving, while the doms' y-component
    #  remains static. Because 'iscontained()' centers the dom, the vals leave
    #  the centered dom much earlier than expected!
    global doms = IntervalBox(domy, tdom.hi..(tdom.hi + tstep))
end

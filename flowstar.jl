using TaylorSeries  # Uses IntervalArithmetic
using TaylorModels  # Version req'd: https://github.com/gaperez64/taylormodels.jl

# Assumes the last variable is t
"""Generate the Taylor polynomial approximation part of the i-th flowpipe
    via Lie derivatives.
"""
function tay_poly(f::Vector{TaylorN{N}}, k::Integer) where {N <: Number}
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
function picard_tm_extension(vector_field_tms::Vector{T}, tmv, domain, k) where T <: TaylorModelN
    # TODO: Does this comment make sense/does it belong here?
    # Represent each component of the actual vector as a TM.
    # The Taylor polynomial approximates the actual function
    # and the remainder represents the error incurred by this
    # approximation.
    f::Vector{TaylorN} = polynomials(vector_field_tms)
    errors::Vector{Interval} = remainders(vector_field_tms)
    # FIXME: Is a copy needed? The substitution f(g(x, t), t)
    # seems to have side effects?
    tmv = copy(tmv)
    # Perform the substitution operation of the Picard operator:
    #   f(g(x, s), s).
    ftm = map((fj) -> evaluate(fj, tmv), f)  # FIXME: Should be a general fun, not TSeries
    # Add the error of aproximating the true vector field.
    add_to_rem(tm, err) = TaylorModelN(polynomial(tm), remainder(tm) + err, tm.x0, domain)
    ftm = add_to_rem.(ftm, errors)
    # Compute the integral's remainder part:
    #   (int_enclosure(pe) + I) * [0, t].
    # FIXME: Once again we assume the t interval is the last one
    tdom = domain.v[end]
    # Prepare two functions to define new error intervals pointwise
    int_enclosure(p_error_terms, dom) = sum(perr(dom) for perr in p_error_terms)
    pe(tm::TaylorModelN) = polynomial(tm)[k:end]
    return map((tmj) ->
	       (int_enclosure(pe(tmj), domain) + remainder(tmj)) * tdom, ftm)
end


"""A broadcasting of the TM remainder function across a TM vector."""
remainders(tmv::Vector{T}) where T <: TaylorModelN = remainder.(tmv)

"""A broadcasting of the TM polynomial function across a TM vector."""
polynomials(tmv::Vector{T}) where T <: TaylorModelN = polynomial.(tmv)

"""Map each component of a box to the degenerate zero interval."""
zerobox(box::IntervalBox)::IntervalBox = IntervalBox(zero(box))

"""Map each component of a box to the degenerate interval of its midpoint."""
midbox(box::IntervalBox)::IntervalBox = IntervalBox(mid(box))

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
function tay_model_error(vector_field_tms::Vector{T},
                         p, domain, k::Integer, J,
                         NR_CONTRACTIVENESS_TRIES::Integer,
                         NR_REFINEMENTS::Integer,
                         REFINEMENT_EPS::Float64,
                         SCALE::Float64;
                         fixpoint_callback::Function=midbox) where T <: TaylorModelN
    # get a copy of the t variable before messing up order
    t = get_variables()[end]

    # To account for composition, we'll up the current order
    old_order = get_order()
    hgr_order = old_order * maximum(get_order.(p))
    vars = set_variables(get_variable_string(), order=hgr_order)

    # Construct the TM vector now.
    # Note: Anonymous function tuple destructuring has unique syntax:
    #   https://discourse.julialang.org/t/argument-destructuring-and-anonymous-functions/24893
    construct_tm = ((poly, rem),) -> TaylorModelN(poly, rem, fixpoint_callback(domain), domain)
    tmv = map(construct_tm, zip(p, J))

    # Add a dummy TM (t, [0, 0]) for the time variable because the substitution
    # requires one TM for every variable, which includes the time variable.
    # FIXME: Can this ugliness be circumvented?
    dummy_t_tm = TaylorModelN(t, 0..0, fixpoint_callback(domain), domain)

    J0 = nothing
    J1 = nothing
    for nct in 1:NR_CONTRACTIVENESS_TRIES
	println("Contractiveness attempt: $nct")
        J0 = remainders(tmv)

        # Perform the first refinement.
        append!(tmv, dummy_t_tm)  # FIXME: appending dummy
        J1 = picard_tm_extension(vector_field_tms, tmv, domain, k)

        # Test contractiveness.
        if all(issubset.(J1, J0))
            # Contractiveness test succeeded, pass along the contractive (safe)
            # remainder J1.
            tmv = map(construct_tm, zip(p, J1))
            break
        # Contractiveness failure condition reached: nr of retries exhausted.
        elseif nct == NR_CONTRACTIVENESS_TRIES
            println("Could not find a contractive remainder.")
            @assert false "Could not find a contractive remainder"
        end

        # Contractiveness test failed, widen all initial remainders.
        tmv = map(construct_tm, zip(p, remainders(tmv) * SCALE))
    end

    for nr in 1:NR_REFINEMENTS
	print("Refinement no. $nr")
        append!(tmv, dummy_t_tm)  # FIXME: appending dummy
        Jn = picard_tm_extension(vector_field_tms, tmv, domain, k)

        max_improvement::Float64 = diam(IntervalBox(remainders(tmv))) - diam(IntervalBox(Jn))
        @assert max_improvement >= 0.0 "Refinement should not decrease the bound tightness!"

        tmv = map(construct_tm, zip(p, Jn))
        
        println("  (max improvement=$max_improvement)")
        if max_improvement < REFINEMENT_EPS
            break
        end
    end

    # Collect the results.
    errors = remainders(tmv)

    # We now restore the order.
    vars = set_variables(get_variable_string(), order=old_order)
    return errors
end

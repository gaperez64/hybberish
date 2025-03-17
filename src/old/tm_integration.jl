using Plots
using TaylorModels
include("flowstar.jl")
include("ode_shifting.jl")

"""An over-approximation of the true flow for a slice of the time horizon.
"""
struct OverApprox
    initial_set::Vector{T} where T <: TaylorModelN
    step_size::Float64
    initial_box::IntervalBox
    updated_box::IntervalBox
end

"""A shorthand for an over-approximation based on the given domains."""
function OverApprox(initial_set::Vector{T},
                    step_size::Float64,
                    dom_initial::IntervalBox,
                    dom_updated::IntervalBox,
                    shift_offsets::Vector{Float64}) where T <: TaylorModelN
    # Applying unshifting here ensures that no user has to be concerned
    # about any shifting that may previously have been involved.
    box_initial = interval_initial_set(initial_set, dom_initial)
    box_initial = unshift_boxes([box_initial], shift_offsets)[1]
    box_initial = IntervalBox(box_initial)
    box_updated = interval_initial_set(initial_set, dom_updated)
    box_updated = unshift_boxes([box_updated], shift_offsets)[1]
    box_updated = IntervalBox(box_updated)
    return OverApprox(initial_set,
                      step_size,
                      box_initial,
                      box_updated)
end

"""Convert the given box to a TM initial set.

    A box also represents a TM initial set. This is a different,
    but equivalent representation to a TM initial set.
"""
function tm_initial_set(box::IntervalBox;
                        fixpoint_callback::Function=midbox)
    # Each variable should have a matching box element.
    @assert(get_numvars() == length(box))

    # Do the conversion: x  =>  (x, [0, 0])  given the domain of x.
    vars = get_variables()
    # FIXME: This also adds an entry for variable t, but it should not?
    #       Because an ODE with components x,y should have a 2-dimensional
    #       init set Xi, with components only for x,y and not t?
    return map((var) -> TaylorModelN(var, 0..0, fixpoint_callback(box), box), vars)
end


"""Convert the given TM initial Xi set to a box Bi.

    A box also represents a TM initial set. This is a different,
    but equivalent representation to a TM initial set.
"""
function interval_initial_set(Xi, domain::IntervalBox)
    # Each variable except t should have a matching TM initial set element.
    @assert(get_numvars() == length(Xi))

    # Do the conversion: (p, I)  =>  p(domains) + I  given the domains.
    Bij(Xij) = polynomial(Xij)(domain) + remainder(Xij)
    return map((Xij) -> Bij(Xij), Xi)
end


"""Given a single TM (p, I), compute the corresponding flowpipe.

    Given an initial set Xi, the flowpipe is computed as (p(Xi, t), I).
"""
function tm_flowpipe(p::TaylorN, I::Interval, Xi::Vector, domain::IntervalBox;
                     fixpoint_callback::Function=midbox)
    Xi = copy(Xi)
    # Compute the flowpipe Fi = (p(Xi), I) = (ps, Is + I)
    #             where p(Xi) = (ps, Is)
    tm = p(Xi)
    ps = polynomial(tm)
    Is = remainder(tm)
    return TaylorModelN(ps, Is + I, fixpoint_callback(domain), domain)
end


"""Compute the next iteration's initial set Xi based on the
    current flowpipe Fi.

    The time step-size δi is used for substitution in a TaylorSeries.
"""
function initial_set(Fi::TaylorModelN, delta_i::TaylorN, variables, domain::IntervalBox;
                     fixpoint_callback::Function=midbox)
    # To compute Xi, fix t=δi in Fi.
    # FIXME: Assume t is the last variable and exclude it from the list.
    vars_except_t = variables[1:end-1]
    valuations = vcat(vars_except_t, [delta_i])
    ps = polynomial(Fi)(valuations)
    Is = remainder(Fi)
    return TaylorModelN(ps, Is, fixpoint_callback(domain), domain)
end

"""Construct a TM vector approximation of the vector field of the ODEs.

	For polynomial ODEs, this results in a vector of TMs with errors 0..0.
	For non-polynomial ODEs, the errors likely differ from 0..0.

	@param[in]        f! The vector field (constructor) function.
	@param[in] variables The variables to use in construction.
	@param[in]   domains The domains to of the variables.
	@return A vector of TMs representing the vector field,
"""
function vector_field_tmv(f!::Function, variables::Vector, domains::IntervalBox;
                          fixpoint_callback::Function=midbox)
	@assert length(variables) == length(domains) "Each variable should specify a domain."

	# FIXME: Suppose t is also part of the variables
	num_vars_no_t::Integer = length(variables) - 1
	# The vector of TMs approximating the vector field components.
	tmv::Vector{TaylorModelN} = Vector{TaylorModelN}(undef, num_vars_no_t)
	# The TM representation of the (shifted) variables.
	# Construct the identity TM for each variable.
	variables_tms::Vector{TaylorModelN} = tm_initial_set(domains, fixpoint_callback=fixpoint_callback)

    # FIXME: Bumping up the order to avoid assertion errors
    # during TM arithmetic.
    old_order = get_order()
    new_order = old_order*2
	_ = set_variables(get_variable_string(), order=new_order)

	# Construct the vector field by callback.
	f!(tmv, variables_tms)

    # TODO: Reset to old order.
    _ = set_variables(get_variable_string(), order=old_order)

	return tmv
end

"""Construct a TM approximation of the vector field of the ODEs.

	For polynomial ODEs, this results in a vector of TMs with errors 0..0.
	For non-polynomial ODEs, the errors likely differ from 0..0.

	@param[in]        f! The vector field (constructor) function.
	@param[in] variables The variables to use in construction.
	@param[in]   domains The domains to of the variables.
	@param[in]     shift If true, then shift the ode's and domains so that
                         the domains are centered on the origin but the overall
						 effect on TM integration is neutral.
                         Optional, the default is false.
	@return (
		A vector of TMs representing the (shifted) vector field,
		A vector of shift offsets with one element for each vector field component,
        The shifted domains
	)
"""
function vector_field_tmv_shifted(f!::Function, variables::Vector, domains::IntervalBox;
                                  shift::Bool=false, fixpoint_callback::Function=midbox)
	@assert length(variables) == length(domains) "Each variable should specify a domain."

	# FIXME: Suppose t is also part of the variables
	num_vars_no_t::Integer = length(variables) - 1
	# The vector of TMs approximating the vector field components.
	tmv::Vector{TaylorModelN} = Vector{TaylorModelN}(undef, num_vars_no_t)
	# The TM representation of the (shifted) variables.
	variables_tms::Vector{TaylorModelN} = Vector{TaylorModelN}(undef, num_vars_no_t)
    # Copy just to be safe, since this may be returned from the function.
    domains_shifted::IntervalBox = IntervalBox(domains)
    shift_offsets::Vector = fill(0.0, length(variables))

	if ! shift
		# Construct the identity TM for each variable.
	    variables_tms = tm_initial_set(domains, fixpoint_callback=fixpoint_callback)
	else
		domains_shifted, variables_tms, shift_offsets =
            shifted_vars_tms(domains, fixpoint_callback=fixpoint_callback)
	end

    # FIXME: Bumping up the order to avoid assertion errors
    # during TM arithmetic.
    old_order = get_order()
    new_order = old_order*2
	_ = set_variables(get_variable_string(), order=new_order)

    f!(tmv, variables_tms)

    # FIXME: Reset to old order.
    _ = set_variables(get_variable_string(), order=old_order)

    return tmv, shift_offsets, domains_shifted
end

"""The standard TM integration algorithm.

    Compute an overapproximation of the true flow of the
    system of ODEs represented by the given vector field,
    as a sequence of flowpipes over partial time horizons
    [0, δi] of the full time horizon [0, Δ].

    @param[in] vector_field_constructor A constructor for the TM vector
                                        representation of the vector field.
    @param[in]                   domain The domains of the ODE variables.
    @param[in]                        k The TM arithmetic/truncation order.
    @param[in]                        J The initial remainder estimate for
                                        finding a contractive remainder.
    @param[in]             time_horizon The bounded time horizon.
    @param[in]           TIME_STEP_SIZE The fixed time step size
    @param[in]       TIME_STEP_SIZE_EPS Any time step size smaller than this
                                        threshold gets skipped. This pertains
                                        to dynamic step sizes, or to handle
                                        edge cases w.r.t. floating point
                                        errors.
    @param[in] NR_CONTRACTIVENESS_TRIES The number of times to attempt widening
                                        the remainder estimate in order to find
                                        a contractive remainder.
    @param[in]           NR_REFINEMENTS The number of additional remainder
                                        refinements to apply after finding the
                                        contractive remainder.
    @param[in]                    SCALE The scale factor by which to widen the
                                        remainder estimate when the
                                        contractiveness test fails.
    @param[in]           REFINEMENT_EPS The minimum improvement a remainder
                                        refinement step should affect. If the
                                        improvement falls below this threshold,
                                        then we declare the refinement to have
                                        converged and stop the refinment loop.
    @param[in]          PRECONDITIONING If true, then apply preconditioning.
    @param[in]   RECOMPUTE_VECTOR_FIELD If true, then reconstruct the TM vector
                                        representation of the ODEs using the
                                        most current variable domains.
    @param[in]        fixpoint_callback A callback to determine the fixpoint
                                        for some Taylor model, given only the
                                        prospective domains of that Taylor
                                        model's variables.
"""
function tm_integration(vector_field_constructor::Function,
                        domain::IntervalBox,
                        k::Integer,
                        J::IntervalBox,
                        time_horizon::Float64,
                        TIME_STEP_SIZE::Float64,
                        TIME_STEP_SIZE_EPS::Float64,
                        NR_CONTRACTIVENESS_TRIES::Integer,
                        NR_REFINEMENTS::Integer,
                        SCALE::Float64;
                        REFINEMENT_EPS::Float64 = 0.001,
                        PRECONDITIONING::Bool = false,
                        RECOMPUTE_VECTOR_FIELD::Bool = false,
                        fixpoint_callback::Function=midbox)
    @assert(time_horizon > 0)  # The time horizon must not be [0, 0].
    @assert(TIME_STEP_SIZE > 0)
    @assert(NR_CONTRACTIVENESS_TRIES >= 0)
    @assert(NR_REFINEMENTS >= 0)

    # Expect the given domains to NOT contain a time component, because
    # the initial hyperrectangle is implicitly situated at t=0.
    domain = IntervalBox(domain..., 0.0..0.0)

    # TODO: Where & how should this comment be integrated into function docs?
    # The vector field f used for TM integration is represented as
    #   f = [ (p1, I1), ... ]
    # a vector of TMs where
    #   pi is a polynomial approximation of the i-th vector field
    #      component, which may have been non-polynomial
    #   Ii is the error accrued by replacing the i-th, true vector
    #      field component with polynomial pi
    # TODO: Currently we're just shifting the ODEs. Properly implement preconditioning!

    vector_field_tms::Vector{T} where T <: TaylorModelN,
    shift_offsets::Vector{Float64},
    domain =
        vector_field_tmv_shifted(vector_field_constructor, get_variables(), domain,
                                 shift=PRECONDITIONING, fixpoint_callback=fixpoint_callback)

    # Derive the corresponding TM initial set from the interval initial set.
    X0 = tm_initial_set(domain, fixpoint_callback=fixpoint_callback)

    # Start the TM integration loop.
    remaining_time::Float64 = time_horizon
    over_approximations::Array{OverApprox} = [
        # Consider the initial conditions an over-approximation.
        OverApprox(X0, 0.0, domain, domain, shift_offsets)
    ]
    dom_i::IntervalBox = domain
    Fi = nothing    # The current flowpipe.
    Xi::Vector = X0 # The current initial set.
    # FIXME: Would a for-loop be cleaner?
    while remaining_time > 0.0
        delta_i = min(remaining_time, TIME_STEP_SIZE)
        remaining_time -= delta_i
        println("time step [0, $delta_i] (remaining Δ: $remaining_time)")

        if delta_i < TIME_STEP_SIZE_EPS
            println("==> Skipping! δi = $delta_i < $TIME_STEP_SIZE_EPS, the step-size is too small.")
            continue
        end

        # Updates the time domain based on the chosen step size?
        # This matters in relation to variable step sizes.
        # FIXME: Assumes the domains contain a time component.
        dom_i = IntervalBox(dom_i) # FIXME: This type casting should not be necessary
        dom_i = IntervalBox(dom_i.v[1:end-1]..., (0.0)..(delta_i))

        if RECOMPUTE_VECTOR_FIELD
            # TODO: Currently we're just shifting the ODEs. Properly implement preconditioning!
            vector_field_tms, shift_offsets, dom_i =
                vector_field_tmv_shifted(vector_field_constructor, get_variables(), dom_i,
                                        shift=PRECONDITIONING, fixpoint_callback=fixpoint_callback)
            # TODO: X_i-1 is used to construct X_i. Note that X_i-1 is already shifted, and we
            #       shift X_i by a different amount. So the shifting applied dom_i is
            #       not enough! The sum of all shifts of all previous iterations should additionally
            #       be applied to dom_i?
            #       Don't forget to pass this changed offset vector the OverApprox constructor as well,
            #       so that unshifting happens correctly!
        end

        println("domains=$dom_i")

        # Step 1: generate the poly approximation of the flow.
        p = tay_poly(polynomials(vector_field_tms), k)
        println("pl=$p")

        # Step2: generate a safe, refined remainder interval.
        I = tay_model_error(vector_field_tms, p, dom_i, k, J,
                            NR_CONTRACTIVENESS_TRIES,
                            NR_REFINEMENTS,
                            REFINEMENT_EPS,
                            SCALE,
                            fixpoint_callback=fixpoint_callback)

        # Step 3: Compute the flowpipe and next initial set.
        # get a copy of the t variable before messing up order
        t = get_variables()[end]
        old_order = get_order()
        old_variables = get_variables()
        # FIXME: AFAIK the 'evaluate' function for TaylorN of TaylorSeries
        # acts up when the substitution values have different orders.
        # For example:
        #   x, y, a, b = set_variables("x y a b", order=3)
        #   p = TaylorN(x + y)
        #   p([1, 2, 3, 4])   # This works fine.
        #   p([a, b, x, y])   # This works fine too.
        #   p([1, b, 3, y])   # ERROR.
        # This even fails if one of the variables is of lower order than others!
        #   x, y, a, b = set_variables("x y a b", order=3)
        #   p = TaylorN(x + y)
        #   p([a, b, x, y])   # This works fine.
        #   a = get_variables(1)[3]
        #   p([a, b, x, y])   # ERROR. a is now order 1, which is less than 3!
        #
        # FIXME: The solution: Convert the scalar δi to a TaylorN.
        # Construct a TaylorN representation of the scalar δi.
        # Must be of the same order as the variables, else evaluation fails!
        delta_i_tayn = TaylorN(delta_i, get_order())

        # To account for composition, we'll up the current order
        hgr_order = old_order * maximum(get_order.(p))
        vars = set_variables(get_variable_string(), order=hgr_order)

        # TODO: Compute the flowpipe components.
        # FIXME: Does TaylorModels actually add the enclosures of truncated
        # terms during TM arithmetic to the remainders? Because the
        # substitution below requires this!
        # Also, does TaylorModels take into account that TaylorSeries seems
        # to compute the order of a constructed expression based on the order
        # of the lowest order term.
        # e.g. for x order 1 and y order 3, then (x + y^2) is order 1 and
        # actually becomes x with y^2 truncated.
        tmv = zip(p, I)
        Fi = map(((pj, Ij),) -> tm_flowpipe(pj, Ij, Xi, dom_i, fixpoint_callback=fixpoint_callback), tmv)

        # TODO: Compute the initial set vector.
        Xi = map((Fij) -> initial_set(Fij, delta_i_tayn, old_variables, dom_i, fixpoint_callback=fixpoint_callback), Fi)

        dummy_t_tm = TaylorModelN(t, 0..0, fixpoint_callback(dom_i), dom_i)
        # Output Xi AFTER appending the dummy TM, for consistency
        # with X0 that was added to the output list before the loop.
        push!(Xi, dummy_t_tm)   # FIXME: Append dummy


        # TODO: This is an exploratory change, does it make sense/work in practice?
        # The new flowpipe Xi represents the current progress of integration.
        # Make the domain the hyperrectangle representation of Xi.
        dom_i = IntervalBox(interval_initial_set(Xi, dom_i))


        # FIXME: Should this OverApprox be created using the dom_i from
        #        BEFORE or AFTER it is updated???
        #    ==> Perhaps AFTER? Since the very first OverApprox already uses
        #        the very first/initial dom_i
        push!(over_approximations, OverApprox(Xi, delta_i, domain, dom_i, shift_offsets))

        # We now restore the order.
        vars = set_variables(get_variable_string(), order=old_order)
    end

    return over_approximations
end

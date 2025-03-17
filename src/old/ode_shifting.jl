using IntervalArithmetic
using TaylorSeries

"""Compute the shifted domains and shifted variables with which to construct the ODEs.

    Given the hyperrectangle B = (b1, ..., bn).
    Then bi is domain of ODE variable xi. If we center the hyperrectangle B
    on the origin, C = (c1, ..., cn) where ci = bi - mid(bi), then we must shift
    the variables in the opposite direction to counteract this. So the
    shifted variables are yi = xi + mid(bi).

    @param[in] dom The domains of all variables, including time t.
    @return (
        C = (c1, ..., cn) where ci is an interval,
        [(y1, [0, 0]), ..., (yn, [0, 0])] where (y1, [0, 0]) is a Taylor model,
        [mid(b1), ..., mid(bn)]
    )
"""
function shifted_vars_tms(dom::IntervalBox; fixpoint_callback::Function=midbox)
	# Center the initial hyper rectangle on the origin.
	# This applies a SUBTRACTION to the domains.
	offsets::Vector{Float64} = map((domj) -> mid(domj), dom)
	# FIXME: Assume the last component is a dummy interval for time t
	offsets[end] = 0
	dom_shifted = dom - IntervalBox(offsets...)	# The SUBTRACTION

	# Construct the identity TM for each variable.
    vars_tms = tm_initial_set(dom_shifted, fixpoint_callback=fixpoint_callback)

	# Shift the ODEs to compensate for the shifted hyper rectangle
	# This applies an ADDITION to the ODEs that compensates the
	# previous SUBTRACTION. The overall effect of shifting the domains
	# one way and the ODEs the opposite way is NEUTRAL.
	# ==> The ODEs are built using the vector of var TMs, so just
	#     modify that vector
	cte_order = get_order() # We MUST use the current order to avoid accidental truncation
	shifted_vars_tay = map((o) -> TaylorN(o, cte_order), offsets)
	vars_tms = vars_tms + shifted_vars_tay # The ADDITION

	return dom_shifted, vars_tms, offsets
end

"""Translate the given boxes by the given vector.

    @return The translated boxes
"""
function unshift_boxes(boxes::Vector{Vector{Interval{Float64}}}, offsets::Vector{Float64})
	# FIXME: Assume the last component corresponds to time t.
	offsets = offsets[1:(end-1)]
	boxes = map((box) -> box[1:(end-1)], boxes)

	# Undo the SUBTRACTION that the shifting applied to the domain.
	return map((box) -> box + offsets, boxes)
end

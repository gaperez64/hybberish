""" Multivariate Taylor model auxiliaries.

    This module heavily draws from the TaylorModels.jl library by
    Louis Benet and David P. Sanders.
"""

"""Fix the order of both Taylor models' polynomial parts to the same order.

    This operation is out-of-place.
    Given input TMs a, b we compute output TMs a', b'.
    To determine a Taylor model's order, see `get_order(::TaylorModelN)`.

    The truncated terms are bounded and added to the remainder of the
    corresponding Taylor model. The bounding uses Horner's rule.

    @pre domain(a) == domain(b)
    @post get_order(a') == get_order(b') == min(get_order(a), get_order(b)) 
    @return The Taylor models (a', b') with their order fixed.
"""
function fixorder(a::TaylorModelN, b::TaylorModelN)
    @assert domain(a) == domain(b)

    # Base case: the orders are already equal.
    get_order(a) == get_order(b) && return a, b

    order = min(get_order(a), get_order(b))
    apol0, bpol0 = polynomial.((a, b))
    apol, bpol = TaylorSeries.fixorder(apol0, bpol0)

    # Compute a bound for the neglected part of the each polynomial
    # FIXME: The TaylorModels library uses the ´centered_dom(a)´ as domain!
    #        Do we simply assume the expansion point is always the origin?
    dom = domain(a)
    bound_a = bound_truncation(apol0, dom, order) + remainder(a)
    bound_b = bound_truncation(bpol0, dom, order) + remainder(b)

    return TaylorModelN(apol, bound_a, domain(a)),
           TaylorModelN(bpol, bound_b, domain(b))
end

"""Compute the interval bound of polynomial terms of degree exceeding order.

    For example, for a polynomial `x + x^2 + x^3 + x^4` and order 2,
    we compute an evaluation for the subexpression `x^3 + x^4` whose
    terms' degrees exceed 2.

    Note that TaylorSeries' evaluation function is used.
    So Horner's rule is applied to compute tighter bounds.

    @param[in] a The Taylor series to extract terms from.
    @param[in] aux The interval domains to use for polynomial evaluation.
    @param[in] order Evaluate all terms of degree exceeding order.
    @return An interval enclosure (bound) of the selected terms.
"""
function bound_truncation(a::TaylorN, aux::IntervalBox, order::Int)
    # If the max polynomial order is smaller than the truncation order,
    # then there is surely nothing to truncate.
    get_order(a) < order && return zero(aux[1])

    # TODO: Can a deepcopy be prevented? Or is it required for calling
    # TaylorSeries evaluation properly? Or is using a conditional ´map´
    # over the HomogenisPolynomials and constructing a new TaylorModel
    # more efficient? Or using a view?
    res = deepcopy(a)
    res[0:order] .= zero(res[0])
    # Call TaylorSeries' evaluation function which implements Horner's rule.
    return res(aux)
end



#
# The Following TaylorModelN auxiliaries are simple wrappers around the TaylorN
# polynomial member's equivalents for these functions.
#
@inline iterate(a::TaylorModelN, state=0) = iterate(polynomial(a), state)
@inline firstindex(a::TaylorModelN) = firstindex(polynomial(a))
@inline lastindex(a::TaylorModelN) = lastindex(polynomial(a))
@inline eachindex(a::TaylorModelN) = eachindex(polynomial(a))
@inline length(a::TaylorModelN) = length(polynomial(a))
@inline size(a::TaylorModelN) = size(polynomial(a))
getindex(a::TaylorModelN, n::Int) = getindex(polynomial(a), n)

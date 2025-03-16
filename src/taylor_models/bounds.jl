""" Multivariate Taylor model bounding functions.

    This module heavily draws from the TaylorModels.jl library by
    Louis Benet and David P. Sanders.
"""

""" Bound the absolute remainder of the polynomial approximation of `f` given
    by the Taylor polynomial `polf` around `x0` on the interval `I`. It requires
    the interval extension `polfI` of the polynomial that approximates `f` for
    the whole interval `I`, in order to compute the Lagrange remainder.

    If `polfI[end]` has a definite sign, then it is monotonic in the intervals
    [I.lo, x0] and [x0.hi, I.hi], which is exploited; otherwise, it is used
    to compute the Lagrange remainder. This corresponds to Prop 2.2.1 in Mioara
    Joldes PhD thesis (pp 52).
"""
function bound_remainder(
        f::Function, polf::Taylor1, polfI::Taylor1,
        x0::Interval{T}, I::Interval{T}) where {T}
    _order = get_order(polf) + 1
    # Determine monotonicity from the final term of the
    # function's (finite order) Taylor interval extension.
    fTIend = polfI[_order]
    bb = sup(fTIend) < 0 || inf(fTIend) > 0
    return _monot_bound_remainder(Val(bb), f, polf, polfI, x0, I)
end

""" Computes the remainder exploiting monotonicity;
    see Prop 2.2.1 in Mioara Joldes' PhD thesis (pp 52).
"""
@inline function _monot_bound_remainder(
        ::Val{true}, f::Function, polf::Taylor1, polfI::Taylor1,
        x0::Interval{T}, I::Interval{T}) where {T}
    # The absolute remainder is monotonic.
    a = Interval(inf(I))
    b = Interval(sup(I))
    rem_lo = f(a) - polf(a-x0)
    rem_hi = f(b) - polf(b-x0)
    rem_x0 = f(x0) - polf[0]
    return hull(rem_lo, rem_x0, rem_hi)
end

""" Computes the remainder using Lagrange bound.
"""
@inline function _monot_bound_remainder(
        ::Val{false}, f::Function, polf::Taylor1, polfI::Taylor1,
        x0::Interval{T}, I::Interval{T}) where {T}
    _order = get_order(polf) + 1
    fTIend = polfI[_order]
    return fTIend * (I-x0)^_order
end

""" Multivariate Taylor model rpa functions.

    This module heavily draws from the TaylorModels.jl library by
    Louis Benet and David P. Sanders.
"""

# Elementary functions
fnlist = (
    :inv,  :sqrt, :exp, :log,
    :sin,  :cos,  :tan,
    :asin, :acos, :atan,
    :sinh, :cosh, :tanh
)

for fn in fnlist
    @eval $fn(tm::TaylorModelN) = rpa($fn, tm)
end

""" Rigorous polynomial approximation (RPA) for the function `g` using the
    Taylor model with absolute/relative remainder `tmf`. The bound is computed
    exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

    The function `f` should support the following call signatures.
        f(::Taylor1{T})
        f(::Taylor1{Interval{T}})
    This interface mimics the elementary functions.
        sin(x)
        cos(x)
        ...
"""
function rpa(g::Function, tmf::TaylorModelN{N,T,S}) where {N,T,S}
    _order::Integer = get_order(tmf)

    # FIXME: Is this commented section relevant?
    # # Avoid overestimations
    # if tmf == TaylorModelN(constant_term(tmf.pol), _order, expansion_point(tmf), domain(tmf))
    #     # ... in case `tmf` is a simple constant polynomial
    #     range_g = (g(tmf.pol))(centered_dom(tmf)) + remainder(tmf)
    #     return TaylorModelN(range_g, _order, expansion_point(tmf), domain(tmf))
    # else
    #     v = get_variables(T, _order)
    #     any( tmf.pol .== v ) && _rpaar(g, expansion_point(tmf), domain(tmf), _order)
    # end

    f_pol0::T = constant_term(tmf)

    # Compute the range of Taylor model tmf, including its remainder, i.e.
    # evaluate an interval enclosure of the Taylor model over its domain.
    # A shorthand for evaluating a TM over its domain is `(::TaylorModelN)()`.
    range_tmf::Interval = tmf()

    # Compute RPA for `g` around tmf's constant term, over range_tmf.
    polg::Taylor1,
    remg::Interval = _rpa(g, f_pol0, range_tmf, _order)

    # Use original independent variable
    tm1 = tmf - f_pol0   # OVER-ESTIMATION; IMPROVE
    tmres = polg( tm1 )

    return TaylorModelN(
        polynomial(tmres),
        remainder(tmres) + remg,
        domain(tmf)
    )
end

""" Rigorous polynomial approximation (RPA) with absolute/relative remainder
    for the function `f` on the interval `I`,  using a Taylor expansion
    around the *interval* `x0` of order `_order`. The bound is computed
    by `bound_remainder`(@ref) exploiting monotonicity if possible, otherwise,
    it uses Lagrange bound.
"""
function _rpa(f::Function, x0::T, I::Interval{T}, _order::Integer) where {T}
    # Taylorise the function via function composition and
    # order k Taylor model arithmetic.
    # = f( x0 +     1*t + O(t^(k+1)) )
    polf::Taylor1  = f( Taylor1([x0, one(x0)], _order  ) )
    # = f(  I + [1,1]*t + O(t^(k+2)) )    
    polfI::Taylor1 = f( Taylor1([I,  one(I) ], _order+1) )
    x0I = Interval(x0)

    rem = bound_remainder(f, polf, polfI, x0I, I)
    # FIXME: The original (TaylorModels.jl) implementation constructs
    #        a TaylorModel1 object. We do not use this struct/type.
    #        Does this negatively affect the result?
    # FIXME: The expansion point is chosen to be x0I, an Interval!
    #        How do we combine this with this implementation
    #        that does not store an expansion point?
    #       ==> Shift the coefficients by the midpoint of x0I?
    #       ==> Does TaylorSeries have a convenience function for shifting those?
    # return TaylorModel1(polf, rem, x0I, I) # FIXME: uncomment original code?

    # We do not implement the TaylorModel1 struct, just return its components.
    return polf, rem
end

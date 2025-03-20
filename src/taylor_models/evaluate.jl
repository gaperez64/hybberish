""" Multivariate Taylor model evaluation in N independent variables.

    This module heavily draws from the TaylorModels.jl library by
    Louis Benet and David P. Sanders.
"""

# TaylorSeries TaylorN evaluate:
#       https://github.com/JuliaDiff/TaylorSeries.jl/blob/0298820a6d1f903185e20849c8178ffb0c1cd503/src/evaluate.jl#L257
# TODO: What is the `sorting` parameter for HomogeneousPolynomial for?
#       https://github.com/JuliaDiff/TaylorSeries.jl/blob/0298820a6d1f903185e20849c8178ffb0c1cd503/src/evaluate.jl#L307
# TODO: Note that `evaluate` != `_evaluate`, mind the underscore.
#       `_evaluate` loops over all HomogeneousPolynomials. Why?
#       https://github.com/JuliaDiff/TaylorSeries.jl/blob/0298820a6d1f903185e20849c8178ffb0c1cd503/src/evaluate.jl#L319
# TODO: When is Horner form / Horner's method applied? See `_horner!`?
#       But `_horner!` is in-place for Taylor1? Maybe it's implemented
#       in-line and out-of-place for TaylorN?



"""Evaluate a TaylorModel on an interval box."""
function evaluate(tm::TaylorModelN{N,T,S}, a::IntervalBox{N,S}) where {N,T,S}
    @assert(iscontained(a, tm),
        "The evaluation values are not subset of the evaluated TaylorModelN's domain."*
        "\na ⊈ domain(tm) : $a ⊈ $(domain(tm))")
    return polynomial(tm)(a) + remainder(tm)
end

"""Evaluate a TaylorModel on an interval box."""
(tm::TaylorModelN{N,T,S})(a::IntervalBox{N,S}) where {N,T,S} = evaluate(tm, a)

"""Broadcast evaluation using the same IntervalBox over a vector of TaylorModelN."""
evaluate(tmv::Vector{TaylorModelN{N,T,S}}, a::IntervalBox{N,S}) where {N,T,S} =
    IntervalBox( [ tmv[i](a) for i in eachindex(tmv) ] )

"""Evaluate a TaylorModel on any vector of values."""
function evaluate(tm::TaylorModelN{N,T,S}, a::AbstractVector{R}) where {N,T,S,R}
    @assert(iscontained(a, tm),
        "The evaluation values are not subset of the evaluated TaylorModelN's domain."*
        "\na ∉ domain(tm) : $a ∉ $(domain(tm))")
    return polynomial(tm)(a) + remainder(tm)
end

"""Evaluate a TaylorModelN on a vector of TaylorModelN values.

   This is a special case of evaluation, it corresponds to composition.
"""
function evaluate(tm::TaylorModelN{N,T,S}, a::Vector{TaylorModelN{N,T,S}})::TaylorModelN{N,T,S} where {N,T,S}
    @assert(iscontained(a, tm),
        "The evaluation values are not subset of the evaluated TaylorModelN's domain."*
        "\na ⊈ domain(tm) : Range(a) ⊈ $(domain(tm))\na = $a")
    composed = polynomial(tm)(a)
    return TaylorModelN(
        polynomial(composed),
        remainder(composed) + remainder(tm),
        domain(composed)
    )
end

"""Evaluate a TaylorModelN using function-call-like syntax."""
(tm::TaylorModelN{N,T,S})(a::AbstractVector{R}) where {N,T,S,R} = evaluate(tm, a)

"""Evaluate a TaylorModelN over its entire domain using function-call-like syntax."""
(tm::TaylorModelN{N,T,S})() where {N,T,S} = evaluate(tm, domain(tm))

"""Substitute a TaylorModelN into a Taylor1.

    Note that the Taylor1 is not a Talor model, and consequently
    it **does not** include a remainder part.

    @param[in] polg The Taylor1 polynomial to substitute into.
    @param[in]  tmf The value to substitute.
    @return The substitution result.
"""
function _evaluate(polg::Taylor1{T}, tmf::TaylorModelN{N,T,S}) where{N,T,S}
    _order = get_order(tmf)
    @assert _order == get_order(polg)

    tmres = TaylorModelN(zero(constant_term(polg)), _order, domain(tmf))
    @inbounds for k = _order:-1:0
        tmres = tmres * tmf
        tmres = tmres + polg[k]
    end

    # Returned result does not include a remainder related to polg.
    return tmres
end

(tm::Taylor1)(x::TaylorModelN) = _evaluate(tm, x)

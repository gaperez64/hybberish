""" Multivariate Taylor model promotion.

    This module heavily draws from the TaylorModels.jl library by
    Louis Benet and David P. Sanders.
"""

function promote(a::TaylorModelN{N,T,S}, b::R) where {N, T<:Real, S<:Real, R<:Real}
    a_prom, b_prom = promote(polynomial(a), b)
    return (
        TaylorModelN(a_prom, remainder(a), domain(a)),
        TaylorModelN(b_prom, zero(remainder(a)), domain(a))
    )
end

promote(b::R, a::TaylorModelN{N,T,S}) where {N, T<:Real, S<:Real, R<:Real} =
    reverse(promote(a,b))

function promote(a::TaylorModelN{N,T,S}, b::TaylorN{R}) where {N, T, S, R}
    RR = TaylorN{promote_type(T, R)}
    return (
        TaylorModelN(convert(RR, polynomial(a)), remainder(a), domain(a)),
        TaylorModelN(convert(RR, b), 0..0, domain(a))
    )
end

promote(b::TaylorN{R}, a::TaylorModelN{N,T,S}) where {N, T, S, R} =
    reverse( promote(a, b) )

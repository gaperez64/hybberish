""" Multivariate Taylor model constructors.

    This module heavily draws from the TaylorModels.jl library by
    Louis Benet and David P. Sanders.
"""

using IntervalArithmetic
using TaylorSeries


""" A Taylor model in N independent variables.

    AbstractSeries is an abstract type defined by the TaylorSeries library.
    It represents a parameterized polynomial type.

    The parameterized types mean the following:
        N    The number of independent variables to consider.
        T    The Taylor polynomial coefficient type.
        S    The remainder interval bounds' type.
"""
struct TaylorModelN{N,T,S} <: AbstractSeries{T}
    pol  :: TaylorN{T}        # polynomial approx (of order `ord`)
    rem  :: Interval{S}       # remainder
    dom  :: IntervalBox{N,S}  # interval of interest

    # Inner constructor
    function TaylorModelN{N,T,S}(
            pol::TaylorN{T},
            rem::Interval{S},
            dom::IntervalBox{N,S}) where {N, T <: TS.NumberNotSeries, S <: Real}

        # The Taylor model requires details about all N independent variables.
        @assert N == get_numvars()
        # If the remainder interval does not contain 0, then the Taylor model's
        # polynomial part falls out of the region specified by `pol + rem`.
        @assert zero(S) in rem

        # FIXME: We currently do not care about the expansion point x0.
        # @assert x0 ⊆ dom

        return new{N,T,S}(pol, rem, dom)
    end
end

#
# Outer constructors
#

"""The basic TaylorModelN constructor."""
TaylorModelN(pol::TaylorN{T}, rem::Interval{S}, dom::IntervalBox{N,S}) where {N,T,S} =
    TaylorModelN{N,T,S}(pol, rem, dom)

# Short-cut for a constant
TaylorModelN(a::Interval{T}, ord::Integer, dom::IntervalBox{N,T}) where {N,T} =
    TaylorModelN(TaylorN(a, ord), zero(dom[1]), dom)
TaylorModelN(a::T, ord::Integer, dom::IntervalBox{N,T}) where {N,T} =
    TaylorModelN(TaylorN(a, ord), zero(dom[1]), dom)

#
# Getter functions.
#

@inline polynomial(tm::TaylorModelN)::TaylorN = tm.pol
@inline remainder(tm::TaylorModelN)::Interval = tm.rem
@inline domain(tm::TaylorModelN)::IntervalBox = tm.dom
@inline get_order(tm::TaylorModelN) = get_order(polynomial(tm))

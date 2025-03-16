""" Multivariate Taylor model arithmetic in N independent variables.

    This module heavily draws from the TaylorModels.jl library by
    Louis Benet and David P. Sanders.
"""

#
# Constants.
#

zero(a::TaylorModelN) =
    TaylorModelN(zero(polynomial(a)), zero(remainder(a)), domain(a))
one(a::TaylorModelN) =
    TaylorModelN(one(polynomial(a)), zero(remainder(a)), domain(a))



#
# Boolean operators.
#

==(a::TaylorModelN, b::TaylorModelN) =
    polynomial(a) == polynomial(b) &&
    remainder(a) == remainder(b) &&
    domain(a) == domain(b)

iscontained(a, tm::TaylorModelN) =
    a in domain(tm)

iscontained(a::IntervalBox{N,S}, tm::TaylorModelN{N,T,S}) where {N,T,S} =
    issubset(a, domain(tm))

iscontained(tmv::Vector{TaylorModelN{N,T,S}}, b::TaylorModelN{N,T,S}) where {N,T,S} =
    iscontained(IntervalBox([ tm(domain(tm)) for tm in tmv ]), b)



#
# Addition and substraction
#

for op in (:+, :-)
    @eval begin
        # Unary + or - operator
        $(op)(a::TaylorModelN) =
            TaylorModelN($(op)(polynomial(a)), $(op)(remainder(a)), domain(a))

        # Binary + or - operator for two Taylor models.
        function $(op)(a::TaylorModelN, b::TaylorModelN)
            @assert domain(a) == domain(b)
            a, b = fixorder(a, b)
            return TaylorModelN(
                $(op)(polynomial(a), polynomial(b)),
                $(op)(remainder(a), remainder(b)),
                domain(a))
        end

        # Bi\nary + or - operator for one Taylor model and a number constant.
        $(op)(a::TaylorModelN, b::T) where {T <: TS.NumberNotSeries} =
            TaylorModelN($(op)(polynomial(a), b), remainder(a), domain(a))
        $(op)(b::T, a::TaylorModelN) where {T <: TS.NumberNotSeries} =
            TaylorModelN($(op)(b, polynomial(a)), $(op)(remainder(a)), domain(a))
    end
end



#
# Multiplication
#

"""The product of two Taylor models."""
function *(a::TaylorModelN, b::TaylorModelN)
    @assert domain(a) == domain(b)
    a_order = get_order(a)
    b_order = get_order(b)
    # The upper bound order for all terms resulting from the product.
    rnegl_order = a_order + b_order
    # TaylorSeries uses a fixed/global, maximum polynomial degree.
    # To make use of the TaylorSeries evaluation function, which
    # implements Horner's rule for polynomial evaluation (!), to bound the
    # truncated terms we must be able to represent the truncated terms as a
    # polynomial object! But the very fact that those terms are of a degree
    # exceeding the maximum order implies that we cannot represent them as
    # such. Thus, to conform to TaylorSeries, and be able to represent the
    # truncated terms, we must limit the input Taylor model orders.
    @assert rnegl_order <= get_order()
    aux = domain(a)


    #
    # Compute the polynomial product.
    #
    a_pol = polynomial(a)
    b_pol = polynomial(b)
    product_pol = a_pol * b_pol # The product is subject to `fixorder(a, b)`
    product_order = get_order(product_pol)
    # The TaylorSeries operator `*(::TaylorN, ::TaylorN)` uses `fixorder`,
    # so the order assigned to the resulting polynomial is the minimum of
    # the two operands' orders. Enforce this post-condition explicitly.
    @assert(product_order == min(get_order(a), get_order(b)))
    # Determine how much `fixorder` lowered the product result's order.
    # = order(a) + order(b) - min(order(a), order(b))
    lowered_orders = rnegl_order - product_order


    #
    # Manually compute an interval bound of all terms that were implicitly
    # dropped due to `fixorder`.
    #
    truncated = Array{HomogeneousPolynomial{TS.numtype(product_pol)}}(
        undef, lowered_orders)
    bound_negl_arr = Array{promote_type(TS.numtype(product_pol), TS.numtype(domain(a)))}(
        undef, lowered_orders)
    # For each truncated order ...
    for truncated_order in product_order+1:rnegl_order
        idx = truncated_order - product_order
        # Initialize the polynomial to the constant 0.
        truncated[idx] = HomogeneousPolynomial(zero(TS.numtype(product_pol)), truncated_order)
        # ... determine all of the product's terms of a given order.
        # = a[0]b[t] + ... + a[i]b[t-i] + ... + a[t]b[0]
        @inbounds for i::Int = 0:truncated_order
            j::Int = truncated_order - i
            # Skip if an out-of-bounds indexation error would occur.
            # Intuitively, an out-of-bounds index corresponds to a
            # HomogeneousPolynomial that is a constant 0 polynomial,
            # so the product would be a constant 0 polynomial as well.
            (i > a_order || j > b_order) && continue
            # Aggregate the polynomial product terms of the given order.
            TaylorSeries.mul!(truncated[idx], a_pol[i], b_pol[j])
        end
        # FIXME: Does the TaylorSeries `evaluate` function for HomogeneousPolynomial use
        # Horner's rule for polynomial evaluation???
        # See:
        #       https://github.com/JuliaDiff/TaylorSeries.jl/blob/0298820a6d1f903185e20849c8178ffb0c1cd503/src/evaluate.jl#L156
        #       https://github.com/JuliaDiff/TaylorSeries.jl/blob/0298820a6d1f903185e20849c8178ffb0c1cd503/src/evaluate.jl#L175
        bound_negl_arr[idx] = truncated[idx](aux)
    end

    # Compute the Taylor model product's remainder part.
    # FIXME: Does the sorting relate to one possible Horner form?
    bound_negl = sum( sort!(bound_negl_arr, by=abs2) )
    bound = remainder_product(a, b, aux, bound_negl)
    product = TaylorModelN(product_pol, bound, domain(a))

    # The resulting TM must propagate the lowest operand polynomial order.
    @assert get_order(product) == product_order
    return product
end

"""Compute the remainder part of the Taylor model product formula.

    @param[in] a The left operand of the product.
    @param[in] b The right operand of the product.
    @param[in] aux The variable domains.
    @param[in] bound_negl The interval enclosure of the truncated (neglected)
                          terms of the product.
    @return The interval remainder part of the Taylor model product formula.
"""
@inline function remainder_product(a::TaylorModelN, b::TaylorModelN,
        aux::IntervalBox, bound_negl::Interval)::Interval
    bound_a = polynomial(a)(aux)
    bound_b = polynomial(b)(aux)
    a_rem = remainder(a)
    b_rem = remainder(b)
    # Apply the Taylor model product remainder formula:
    # = Int(p1) · I2 + I1 · Int(p2) + I1 · I2 + Int(pe)
    return bound_a * b_rem + bound_b * a_rem + a_rem * b_rem + bound_negl
end

"""The product of a Taylor model with a scalar number."""
*(b::T, a::TaylorModelN) where {T <: TS.NumberNotSeries} =
    TaylorModelN(polynomial(a) * b, b * remainder(a), domain(a))

"""The product of a Taylor model with a scalar number."""
*(a::TaylorModelN, b::T) where {T <: TS.NumberNotSeries} = b * a



#
# Division
#

const ERR_TMN_DIV = "Unimplemented; TaylorModelN division is redundant as "*
                    "we Taylorize the dynamics which removes all divisions."

/(a::TaylorModelN, b::TaylorModelN) = throw(ERR_TMN_DIV)
/(b::T, a::TaylorModelN) where {T <: TS.NumberNotSeries} = throw(ERR_TMN_DIV)
/(a::TaylorModelN, b::T) where {T <: TS.NumberNotSeries} = throw(ERR_TMN_DIV)



#
# Powers
#

""" Make explicit that TaylorModelN exp (^) is handled by Base functions
    by mimicking the Base exp (^) function.

    We depend on the following Base library functions.
    - `Base.literal_pow(f::typeof(^), x, ::Val{p})`
    - `Base.:^(x::Number, p::Integer)`
    - `Base.power_by_squaring(x_, p::Integer; mul=*)`

    Since we overload `Base.:^(::Number, ::Integer)` with an implementation
    for a subtype, i.e. `TaylorModelN <: Number`, julia's dispatch mechanism
    does not allow us to call the overloaded function directly. We must instead
    mimic the overloaded function's implementation, but for the subtype.
"""
function ^(x::TaylorModelN, p::Integer)
    # Assert the supertype to show that we mimic the Base exp (^) function:
    # `Base.:^(x::Number, p::Integer)`
    @assert typeof(x) <: Number
    # Make explicit the dependence on the Base library function.
    return Base.power_by_squaring(x, p)
end

"""This file contains examples for our adaption of the TaylorModels.jl library.
"""

include("../taylor_models/BasicTaylorModels.jl")

#
# Test TaylorModelN construction.
#
const ORDER = 3
const REMAINDER = -0.5..0.5
const DOMAIN = IntervalBox(-1..1, -2..2)
x, y = set_variables("x y", order=ORDER)
const POLYNOMIAL = x + y + x*y
tn = TaylorModelN(POLYNOMIAL, REMAINDER, DOMAIN)



#
# Test TalorModelN getters.
#
@assert(ORDER == get_order(tn))
@assert(POLYNOMIAL == polynomial(tn))
@assert(REMAINDER == remainder(tn))
@assert(DOMAIN == domain(tn))
println("## Test TaylorModelN getters")
println("get_order(tn): $(get_order(tn))")
println("polynomial(tn): $(polynomial(tn))")
println("remainder(tn): $(remainder(tn))")
println("domain(tn): $(domain(tn))")
println()



#
# Test TaylorModelN ´fixorder´.
#

# Every Taylor model displays their order via Big-O representation!
displayBigO(true)
println("## Test TaylorModelN fixorder")

x, y = set_variables("x y", order=3)
tx::TaylorModelN = TaylorModelN(x, 0..0, IntervalBox(-1..1, -2..2))
ty::TaylorModelN = TaylorModelN(y, 0..0, IntervalBox(-1..1, -2..2))


# If both TaylorModelN objects have the same order, they remain unchanged.
@assert(get_order(tx) == get_order(ty))
txres, tyres = fixorder(tx, ty)
@assert(tx == txres && ty == tyres)
println("fixorder($tx, $ty)")
println("      = ($txres, $tyres)")

# The order of the objects differs, lower both of them to be the same order.
# This may result in truncation + bounding of the truncated part.
# Define two polynomials:
#   a = x     of order 2
#   b = y^3   of order 3
xlow = get_variables(2)[1] # Lower x's order
ypow3 = y^3
tx = TaylorModelN(xlow,  0..0, IntervalBox(-1..1, -2..2))
ty = TaylorModelN(ypow3, 0..0, IntervalBox(-1..1, -2..2))

# We define the objects to have different order.
@assert(get_order(xlow) != get_order(ypow3))
txres, tyres = fixorder(tx, ty)
# Now the orders should be equalized.
@assert(get_order(txres) == get_order(tyres))
# Lowering the order of polynomial b to 2 truncates its
# singular term "y^3", so b = 0 because no terms remain.
@assert(xlow == polynomial(txres) && 0.0 == polynomial(tyres))
# No truncation means no change to the remainder.
@assert(remainder(tx) == remainder(txres))
# The interval enclosure of the truncated terms is added to the remainder.
tgt_rem = remainder(ty) + polynomial(ty)(domain(ty))
@assert(tgt_rem == remainder(tyres))

println("fixorder($tx, $ty)")
println("      = ($txres, $tyres)")
println()



#
# Test TaylorModelN auxiliary functions:
# indexation, iteration and printing, etc.
#

# Test TaylorModelN `show` and `pretty_print`.
println("## Test TaylorModelN printing")
println(tn)
println()

# Test TaylorModelN `iterate`.
println("## Test TaylorModelN iteration")
for e in tn
    println(e)
end
println()

# Test TaylorModelN `getindex`.
println("## Test TaylorModelN index access")
for idx in eachindex(tn)
    @assert(tn[idx] == polynomial(tn)[idx])
    println("getindex(tn, $idx): $(tn[idx])")
end
println()

# Test remaining auxiliary functions.
println("## Test TaylorModelN axiliaries")
@assert(firstindex(tn) == 0)
@assert(lastindex(tn) == ORDER)
@assert(eachindex(tn) == 0:ORDER)
@assert(length(tn) == ORDER+1)
@assert(size(tn) == (ORDER+1,))
println("firstindex(tn): $(firstindex(tn))")
println("lastindex(tn): $(lastindex(tn))")
println("eachindex(tn): $(eachindex(tn))")
println("length(tn): $(length(tn))")
println("size(tn): $(size(tn))")
println()



#
# Specify which exponentiation implementation is used:
# a custom `^` function or `Base.#^`?
#

# Setup.
println("## Check if there is custom TaylorModelN exponentiation (^)")

xx, yy = set_variables("x y", order=16)
xx, yy = get_variables(8)
tx = TaylorModelN(xx, -0.5..0.5, IntervalBox(-0.5..0.5, -4..4))
ty = TaylorModelN(xx*yy, -1..1,  IntervalBox(-0.5..0.5, -4..4))


#= Note the following call stack of functions in `Base.intfuncs.jl`.

  (1)   Base.literal_pow(f::typeof(^), x, ::Val{p})
  (2)   Base.:^(x::Number, p::Integer)
  (3)   Base.power_by_squaring(x_, p::Integer; mul=*)

  Function (1) represents expression x^p where p is a strictly positive
  literal exponent and x is the Taylor model. Note that the chosen literal_pow
  implementation does not have type restrictions on x, which is why this
  version is called on the TaylorModelN x; it is the fallthrough implementation
  of literal_pow.

  Function (2) requires `x <: Number`. Note that a TaylorModelN object x
  satisfies `x <: AbstractSeries{T <: Number} <: Number` by definition of
  the TaylorModelN struct.

  Function (3) implements the literal exponentiation procedure by unrolling it
  into a sequence of multiplications. By defining the Taylor model product `*`
  we avoid the need to also define the exponent operator `^`.
=#

# Specify some pre-consitions / assumptions.
@assert(typeof(tx) <: Number) # TaylorModelN <: Number

# We assume a custom TaylorModelN exponentiation DOES NOT exist; we expect
# a call stack of Base functions to implement exponentiation.
result = tx^4
@assert(result == Base.literal_pow(^, tx, Base.Val(4)))
@assert(result == Base.:^(tx, 4))
@assert(result == Base.power_by_squaring(tx, 4))

# Test a more complex expression, where a function encodes an ODE expression.
"""A simple ode function."""
ode(x, y) = 4 + x^4 * y
result = ode(tx, ty)
actual = 4 + Base.literal_pow(^, tx, Base.Val(4)) * ty
@assert(result == actual)

println()



#
# Test TaylorModelN arithmetic.
#

# Setup.
set_variables("x y", order=ORDER)
tn = TaylorModelN(POLYNOMIAL, REMAINDER, DOMAIN)

# Test unary + operator.
println("## Test unary + operator.")
res = +tn
@assert(polynomial(tn) == polynomial(res))
@assert(remainder(tn) == remainder(res))
@assert(domain(tn) == domain(res))
println("
  + $tn
= $res
")

# Test binary + operator.
println("## Test binary + operator: TaylorModelN + TaylorModelN.")
res = tn + tn
tgt_pol = 2*x + 2*y + 2*x*y
tgt_rem = remainder(tn) + remainder(tn)
@assert(tgt_pol == polynomial(res))
@assert(tgt_rem == remainder(res))
@assert(domain(tn) == domain(res))
println("
  $tn +
  $tn
= $res
")

# Test binary + operator.
cte = 4
println("## Test binary + operator: TaylorModelN + $cte.")
res = tn + cte
tgt_pol = 4 + x + y + x*y
@assert(tgt_pol == polynomial(res))
@assert(remainder(tn) == remainder(res))
@assert(domain(tn) == domain(res))
println("
  $tn +
  $cte
= $res
")

# Test binary + operator.
cte = 4
println("## Test binary + operator: $cte + TaylorModelN.")
res = cte + tn
tgt_pol = 4 + x + y + x*y
@assert(tgt_pol == polynomial(res))
@assert(remainder(tn) == remainder(res))
@assert(domain(tn) == domain(res))
println("
  $cte +
  $tn
= $res
")

# Test unary - operator.
println("## Test unary - operator.")
tgt_pol = -x - y - x*y
tgt_rem = - remainder(tn)
res = -tn
@assert(tgt_pol == polynomial(res))
@assert(tgt_rem == remainder(res))
@assert(domain(tn) == domain(res))
println("
  - $tn
= $res
")

# Test binary - operator.
println("## Test binary - operator: TaylorModelN - TaylorModelN.")
res = tn - tn
tgt_pol = 0
tgt_rem = remainder(tn) - remainder(tn)
@assert(tgt_pol == polynomial(res))
@assert(tgt_rem == remainder(res))
@assert(domain(tn) == domain(res))
println("
  $tn -
  $tn
= $res
")

# Test binary - operator.
cte = 4
println("## Test binary - operator: TaylorModelN - $cte.")
res = tn - cte
tgt_pol = -4 + x + y + x*y
@assert(tgt_pol == polynomial(res))
@assert(remainder(tn) == remainder(res))
@assert(domain(tn) == domain(res))
println("
  $tn - $cte
= $res
")

# Test binary - operator.
cte = 4
println("## Test binary - operator: $cte - TaylorModelN.")
res = cte - tn
tgt_pol = 4 - x - y - x*y
tgt_rem = - remainder(tn)
@assert(tgt_pol == polynomial(res))
@assert(remainder(tn) == remainder(res))
@assert(domain(tn) == domain(res))
println("
  $cte -
  $tn
= $res
")

println()



#
# Test fixorder in TaylorModelN arithmetic.
#

# Every Taylor model displays their order via Big-O representation!
displayBigO(true)
set_variables("x y", order=ORDER)
tn = TaylorModelN(POLYNOMIAL, REMAINDER, DOMAIN)
println("## Test TaylorModelN fixorder in arithmetic")

# The order of the objects differs, lower both of them to be the same order.
# This may result in truncation + bounding of the truncated part.
# Define two polynomials:
#   a = x     of order 2
#   b = y^3   of order 3
xpol = get_variables(2)[1] # Lower x's order
ypol = y^3
tx = TaylorModelN(xpol, 0..0, IntervalBox(-1..1, -2..2))
ty = TaylorModelN(ypol, 0..0, IntervalBox(-1..1, -2..2))
cte = 4


# Test fixorder in binary + operator.
println("## Test binary + operator: TaylorModelN + TaylorModelN.")
res = tx + ty
tgt_pol = polynomial(tx)
tgt_rem = remainder(tx) + (remainder(ty) + polynomial(ty)(domain(ty)))
# The order is pulled down to that lowest one.
@assert(get_order(tx) < get_order(ty))
@assert(get_order(res) == min(get_order(tx), get_order(ty)))
@assert(tgt_pol == polynomial(res))
@assert(tgt_rem == remainder(res))
@assert(domain(tx) == domain(res))
println("
  $tx +
  $ty
= $res
")

# Test fixorder in binary + operator.
# Constants only contribute to the polynomial part.
println("## Test binary + operator: TaylorModelN + $cte.")
res = tx + cte
tgt_pol = polynomial(tx) + cte
# The constant does not affect the order.
@assert(get_order(res) == get_order(tx))
@assert(tgt_pol == polynomial(res))
@assert(remainder(tx) == remainder(res))
@assert(domain(tx) == domain(res))
println("
  $tx +
  $cte
= $res
")

# Test fixorder in binary - operator.
println("## Test binary - operator: TaylorModelN - TaylorModelN.")
res = tx - ty
tgt_pol = polynomial(tx)
tgt_rem = remainder(tx) - (remainder(ty) - polynomial(ty)(domain(ty)))
# The order is pulled down to that lowest one.
@assert(get_order(tx) < get_order(ty))
@assert(get_order(res) == min(get_order(tx), get_order(ty)))
@assert(tgt_pol == polynomial(res))
@assert(tgt_rem == remainder(res))
@assert(domain(tx) == domain(res))
println("
  $tx -
  $ty
= $res
")

# Test fixorder in binary - operator.
# Constants only contribute to the polynomial part.
println("## Test binary - operator: $cte - TaylorModelN.")
res = cte - tx
tgt_pol = cte - polynomial(tx)
tgt_rem = - remainder(tx)
# The constant does not affect the order.
@assert(get_order(res) == get_order(tx))
@assert(tgt_pol == polynomial(res))
@assert(tgt_rem == remainder(res))
@assert(domain(tx) == domain(res))
println("
  $cte -
  $tx
= $res
")

# Test binary * operator.
println("## Test binary * operator: TaylorModelN * TaylorModelN.")
# Setup so that product does not violate its pre-conditions / assertions.
xm = get_variables(1)[1] # order 1
ym = get_variables(2)[2] # order 2
dom = IntervalBox(-1..1, -2..2)
tx = TaylorModelN(1 + xm,        -2..2, dom)
ty = TaylorModelN(1 + ym + ym^2, -4..4, dom)
# Perform the product.
res = tx * ty
#= Write out the polynomial expression:
  polynomial(tx) * polynomial(ty)
  = (1 + x) * (1 + y + y^2)
  = 1 + x + y + x*y + y^2 + x*y^2
  = (1 + x + y)  + (x*y + y^2 + x*y^2)
  pe = x*y + y^2 + x*y^2
  p1*p2 - pe = 1 + x + y
=#
tgt_pol = 1 + xm + ym
#= Write out the interval remainder expression:
  Int(p1) · I2 + I1 · Int(p2) + I1 · I2 + Int(pe)
=#
xi, yi = domain(tx)
tgt_rem =
  polynomial(tx)(domain(tx)) * remainder(ty) +
  remainder(tx) * polynomial(ty)(domain(ty)) +
  remainder(tx) * remainder(ty) +
  (xi*yi + yi^2 + xi*yi^2)
@assert(get_order(res) == min(get_order(tx), get_order(ty)))
@assert(tgt_pol == polynomial(res))
@assert(tgt_rem == remainder(res))
@assert(domain(tx) == domain(ty) == domain(res))
println("
  $tx *
  $ty
= $res
")



#
# Test the time-cost of high polynomial orders.
#

"""Test how long a Taylor series multiplication of a given order takes.

  Test the following expression:

    x^xpow * y^ypow

  for Taylor series arithmetic of order `order`.
  Perform this test `iterations` times and average the results.
"""
function test_xorder_arithmetic(xpow::Int, ypow::Int, order::Int, iterations::Int)
  @assert iterations > 0
  # The `set_variables` call may take long, so it is considered separately.
  # Calling `set_variables` a second time with the same order as the current
  # max order avoids reconstructing the TaylorSeries hash tables, so is faster.
  pre = time()
  x, y = set_variables("x y", order = order)
  post1 = time()

  summ = 0
  post2 = post1
  for _ = 1:iterations
    res = x^xpow * y^ypow

    temp = time()
    summ += temp - post2
    post2 = temp
  end
  println("
 order(x) = $(get_order(x))
 order(y) = $(get_order(y))
max order = $(get_order())
$x * $y = $res
$(round(post1 - pre, digits=8)) seconds to run `set_variables`.
$(round(summ / iterations, digits=8)) seconds avg across $iterations iterations to run the product op.
$(round(post2 - pre, digits=8)) seconds total time taken.
")
end

println("## Test the time-cost of high polynomial orders.")

# Test the a low polynomial order case as a baseline reference.
test_xorder_arithmetic(1, 1, 10, 20)

# Test high a polynomial order with only low order terms:
#    x^1 + O(||x||^301)
#    y^1 + O(||x||^301)
test_xorder_arithmetic(1, 1, 300, 20)

# Test high a polynomial order with high order terms:
#    x^100 + O(||x||^501)
#    y^200 + O(||x||^501)
test_xorder_arithmetic(100, 200, 500, 20)



#
# Test TaylorModelN evaluation.
#

x, y = set_variables("x y", order=4)


# Test eval using a vector of Int as values.
# Trigger ´evaluate(::TaylorModelN{N,T,S}, ::AbstractVector{R})´.
t = TaylorModelN(x + y, -1..1, IntervalBox(-2..2, -3..3))
values = [-1, 2]
expect = interval(values[1]) + interval(values[2]) + remainder(t)
@assert(evaluate(t, values) == expect)
@assert(t(values) == expect) # The functions call interface is equivalent.

# Test eval using an IntervalBox as values
# Trigger ´evaluate(::TaylorModelN{N,T,S}, ::IntervalBox{N,S})´.
t = TaylorModelN(x + y, -1..1, IntervalBox(-2..2, -3..3))
values = IntervalBox(-1..2, 2..3)
expect = values[1] + values[2] + remainder(t)
@assert(evaluate(t, values) == expect)
@assert(t(values) == expect) # The functions call interface is equivalent.

# Test broadcasting eval using an IntervalBox as values over a vector of
# TaylorModelN objects.
# Trigger ´evaluate(::Vector{TaylorModelN{N,T,S}}, ::IntervalBox{N,S})´.
t = TaylorModelN(x + y, -1..1, IntervalBox(-2..2, -3..3))
tv = [t, t, t]
values = IntervalBox(-1..2, 2..3)
expect = values[1] + values[2] + remainder(t)
@assert(all(r -> r == expect, evaluate(tv, values)))

# Test eval using no values, which defaults to the Taylor model domain
# as values.
# Trigger ´evaluate(::TaylorModelN{N,T,S})´.
t = TaylorModelN(x + y, -1..1, IntervalBox(-2..2, -3..3))
values = domain(t)
expect = values[1] + values[2] + remainder(t)
@assert(t() == expect)

# Test eval using a vector of TaylorModelN objects as values.
# Trigger ´evaluate(::TaylorModelN{N,T,S}, ::AbstractVector{R})´.
# for R being ´Vector{TaylorModelN{N,T,S}}´
x4, y4 = set_variables("x y", order=4)
x2, y2 = get_variables(2)
t = TaylorModelN(x4 + y4, -1..1, IntervalBox(-2..2, -3..3))

# The domains of the substitution values must be the same!
# Taylor model arithmetic expects this.
try
    global tx = TaylorModelN(x2, -0.5..0.5, IntervalBox(-0.5..0.5, -2..2))
    global ty = TaylorModelN(x2*y2, -1..1, IntervalBox(-0.5..0.5, -1..1))
    @assert(domain(tx) != domain(ty))
    t([tx, ty])
    @assert(false)
catch
end

# Use the same domains for the substitution values.
# But, we must still watch out for other assertions in Taylor model arithmetic!
try
    # Violate the Taylor model product pre-condition that requires
    # limiting the operand Taylor model orders w.r.t. ´get_order()´.
    @assert(get_order() == 4) # 3+3 = 6 > 4, which will trigger the assertion.
    x3, y3 = get_variables(3)

    global tx = TaylorModelN(x3, -0.5..0.5, IntervalBox(-0.5..0.5, -2..2))
    global ty = TaylorModelN(x3*y3, -1..1, IntervalBox(-0.5..0.5, -1..1))
    @assert(domain(tx) == domain(ty))
    t([tx, ty])
    @assert(false)
catch
end


# Use the same domains for the substitution values.
# This evaluation should actually succeed.
tx = TaylorModelN(x2, -0.5..0.5, IntervalBox(-0.5..0.5, -4..4))
ty = TaylorModelN(x2*y2, -1..1,  IntervalBox(-0.5..0.5, -4..4))
tz = t([tx, ty])
values = [0, 4]
ix = interval(values[1])
iy = interval(values[2])
expect = (-1..1) + (1..1)*ix + (1..1)*ix*iy + remainder(tz)
# The Taylor model resulting from the substitution inherits the common domain
# of the substitution values, not the domain of the Taylor model that is
# substituted into!
@assert(! iscontained(values, t)) # The result does NOT inherit this domain.
@assert(iscontained(values, tz)) # The result inherits the subst. value domain.
@assert(tz(values) == expect)



#
# Test membership, subset, ...
#

x, y = set_variables("x y", order=4)
t = TaylorModelN(x + y, -1..1, IntervalBox(-2..2, -3..3))
tx = TaylorModelN(x, -0.5..0.5, IntervalBox(-0.5..0.5, 0..0))
ty = TaylorModelN(x*y, -1..1, IntervalBox(0.5..0.5, -2..2))

# Trigger ´iscontained(::Any, ::TaylorModelN)´.
@assert(iscontained([0, 0], t))
@assert(iscontained([-2, 3], t))
@assert(! iscontained([-100, 0], t))

# Trigger ´iscontained(::IntervalBox, ::TaylorModelN)´.
@assert(iscontained(IntervalBox(-1..1, -2..2), t))
@assert(iscontained(IntervalBox(-2..2, -3..3), t))
@assert(! iscontained(IntervalBox(-2..2, -3..100), t))

# Trigger ´iscontained(::::Vector{TaylorModelN{N,T,S}}, ::TaylorModelN{N,T,S})´.
# @assert(iscontained(t, t)) # FIXME: Complete this / fix this



#
# Show how to combine TaylorModelN evaluation and TaylorSeries set_variables.
#

const MAX_ORDER = 8
const ACT_ORDER = Int(MAX_ORDER // 2)

# Set the maximum order to 2K.
# This is stored internally by TaylorSeries and is accessible via get_order().
# This is used in assertions by TaylorSeries.
x, y = set_variables("x y", order=MAX_ORDER)

# Actually make use of polynomial order / truncation order 2K/2 = K.
# This compensates for the Taylor model arithmetic product pre-condition, which
# requires:
#     order(TM1) + order(TM2) <= get_order() = 2K
x, y = get_variables(ACT_ORDER)

# FIXME: Complete this / fix this


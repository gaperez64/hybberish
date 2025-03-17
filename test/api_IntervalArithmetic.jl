using IntervalArithmetic


# Intervals are of the form [a, b] = {x in R: a <= x <= b} where R is the
# set of reals.  Invalid intervals, e.g. where a > b as in [3, 1], are mapped
# to the empty set / empty interval.
invalid = 3..1 # Invalid interval definition: a <= b must hold!
@assert(invalid == IntervalArithmetic.emptyinterval())


# An exception is made for (semi) infinite intervals!
# Neither +inf nor -inf are element of [-inf, +inf]!
infinity = -Inf..Inf
pinfin = Inf in infinity
ninfin = (-Inf) in infinity
@assert(! pinfin)
@assert(! ninfin)
println(" $(Inf) in $infinity : $pinfin")
println("$(-Inf) in $infinity : $ninfin")


# Sometimes loose evaluation is applied during interval arithmetic.
# This means discarding part of the interval that is not included in the domain
# of a function. e.g. sqrt is not defined for negative numbers.
int = -2..4
res = sqrt(int)
@assert(res == 0..2) # [-2, 0[ is discarded
println("sqrt($int) = $res because [-2, 0[ was discarded from $int")


# An Interval contains a decoration specifying some properties of the interval.
# The decoration may change by applying a function to it.
div(x) = 1 / x
int = DecoratedInterval(-2, 4)
println("dec of $int: $(decoration(int))")        # com::Decoration = 4
println("dec of 1/$int: $(decoration(div(int)))") # dac::Decoration = 3
# The information loss incurred by pruning [-2, 0[ from the domain [-2, 4]
# results in decoration "trv" or trivial. The pruning was done due to
# "loose evaluation"; [-2, 0[ is not part of the sqrt function's domain.
@assert(decoration(sqrt(int)) == trv)
println("dec of sqrt($int): $(decoration(sqrt(int)))")
# Ill-formed intervals have a corresponding decoration.
# Interval [a, b] violates a <= b.
int = DecoratedInterval(3, 1)
@assert(decoration(int) == ill)
println("dec of $int: $(decoration(int))")
# NaN is not a well-formed interval.
int = DecoratedInterval(NaN)
@assert(decoration(int) == ill)
println("dec of $int: $(decoration(int))")


# FIXME: Bump up the version of IntervalArithmetic until
# `isguaranteed()` is defined. Then uncomment these tests.
if false
    # The `interval()` constructor results in a guaranteed interval.
    @assert(isguaranteed(interval(1.)) == true)
    # Conversion from a float constant to an interval does not
    # result in a guaranteed interval.
    @assert(isguaranteed(convert(Interval{Float64}, 1.)) == false)
    # Conversion from an interval to an interval results in a guaranteed interval.
    @assert(isguaranteed(convert(Interval{Float32}, interval(-1, 1))) == true)
end


# For interval i = [a, b], `mid(i)` computes the interval's midpoint (a+b)/2.
# `mid` returns a number for an Interval.
const deg1::Interval =  1..1
const int1::Interval = -1..1
const int2::Interval = 0.5..1.5
const box::IntervalBox = IntervalBox(deg1, int1, int2)
m = mid(int1)
@assert(m == 0)
println("mid($int1) = $m")
# `mid` returns a static vector of numbers for an IntervalBox.
m = mid(box)
@assert(typeof(m) <: IntervalArithmetic.StaticArraysCore.SVector)
@assert(m == [1.0, 0.0, 1.0])
println("mid($box) = $m")

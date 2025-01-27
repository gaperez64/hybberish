using TaylorSeries  # Uses IntervalArithmetic
using TaylorModels  # Version req'd: https://github.com/gaperez64/taylormodels.jl


# NOTE: Assumes the last variable is t
"""Generate the Taylor polynomial approximation part of a Taylor model for the
   given function and up to the given degree via Lie derivatives.
"""
function tay_poly(f::Vector{TaylorN{N}}, k::Integer) where {N <: Number}
    vars = get_variables()
    t = vars[end]
    # Let's pad f with a 1 at the end for t
    fp1 = copy(f)
    push!(fp1, 1)
    # Also, prepare a valuation vector with t=0
    # NOTE: We're cheating to obtain a zero with same order as the other
    # variables
    val0 = copy(vars)
    val0[end] = t - t
    # Prepare a first lie derivative and the result
    g = copy(vars)
    deleteat!(g, length(vars))
    # Start a vector function for the result
    res = copy(g)
    
    for i = 1:k
        g = TaylorSeries.jacobian(g, vars) * fp1
        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
        res += term
    end
    
    return res
end


# Dynamics
g(y,t) = -y - sin(t) + cos(t)

# Integration task specification
# a. Take delta_t = 1
# b. Construct a flowpipe consisting of 4 Taylor models, including the initial
#    one
# c. Work with order 4
ord = 4
domy = -2..2

# Initial state bounds
# y(0) = [1, 1]
# t(0) = [0, 0]
vals = IntervalBox(1..1, 0..0)

# Taylor variables (from TaylorSeries library)
y, t = set_variables("y t", order=ord)

# Iteration 1
#
# Step 0 - Taylorize the dynamics:
# We want to have a polynomial approximation of the dynamics centered around
# the midpoint of the current values.
ytm = TaylorModelN(y, interval(0), IntervalBox(mid(vals)), IntervalBox(domy, 0..1))
ttm = TaylorModelN(t, interval(0), IntervalBox(mid(vals)), IntervalBox(domy, 0..1))
# FIXME: Hack to allow for higher degree terms in the intermediate computation
y, t = set_variables("y t", order=ord*2)
gtm = g(ytm, ttm)
# FIXME: We go back to the lower degree afterwards
y, t = set_variables("y t", order=ord)
# We will be needing the highest-degree terms of it too
htm = TaylorModelN(TaylorN(polynomial(gtm)[ord], ord),
		   interval(0), expansion_point(gtm), domain(gtm))

println("Taylor model overapproximation of the dynamics:")
println(gtm)
println("Highest degree terms of overapprox:")
println(htm)
println("Range bound of the latter:")
println(quadratic_fast_bounder(htm))

# Step 1 - Obtain the polynomial part of the Taylor model
p = tay_poly([polynomial(gtm) - polynomial(htm)], 4)
println(p)

# Step 2: Obtain the remainder/error interval of the TM
# TODO

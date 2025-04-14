include("../src/tm_integration.jl")


#
# Test the Lie derivative routine AND picard operator routine.
#
# Base this test off a combination of the the naive TM integration example
# and the QR preconditioned TM integration example computations of the
# Taylor models paper by M. Neher et al. from 2006.

k = 3
set_variables("a b t", order=2*k)
vars = get_variables(k)
a, b, t = vars

# These are the normalized variables, so that the space variable domains would
# be [-1, 1]^n at the start of integration.
#       u(0) =  1 + a
#       v(0) = -1 + b
vars_normalized = [ 1 + a, -1 + b ]

#################################
### 1st INTEGRATION ITERATION ###
#################################
# The Taylorized vector field.
# Due to domain normalization, substitute variables u, v with a, b.
#       u' = v   => b
#       v' = u^2 => a^2
fpoly = [ b, a^2 ]

# The initial variable valuation.
# In the 1st iteration, this is computed BEFORE integration started, it is
# the polynomial that represents the zero-th Taylor model initial set.
val0 = vars_normalized

# Compute the Taylor polynomial approx. of the true flow.
p_lie = tay_poly(fpoly, k, vars, val0)
p_pic = picard(fpoly, k, vars, val0)

actual = [
    1 + a - t  + b*t + (1/2)*t^2 + a*t^2 - (1/3)*t^3,
    -1 + b + t + 2*a*t - t^2 + (a^2)*t - a*t^2 + b*t^2 + (2/3)*t^3
]
# The initial valuations do NOT contain floating point numbers, so
# exact equality (==) is even possible against the expected results.
@assert p_lie == actual
@assert p_pic == actual

#################################
### 2nd INTEGRATION ITERATION ###
#################################
# The initial variable valuation.
# In the 2nd iteration, this is the the previous iteration's integration result.
val0 = [
    0.904667 + 1.01*a + 0.1*b,
    -0.909333 + 0.19*a + 1.01*b + 0.1*(a^2)
]

# Compute the Taylor polynomial approx. of the true flow.
p_lie = tay_poly(fpoly, k, vars, val0)
p_pic = picard(fpoly, k, vars, val0)

actual = [
    0.904667 + 1.01*a + 0.1*b - 0.909333*t + 0.19*a*t + 1.01*b*t +
        0.409211*(t^2) + 0.1*(a^2)*t + 0.913713*a*(t^2) + 0.0904667*b*(t^2) -
        0.274215*(t^3),
    -0.909333 + 0.19*a + 1.01*b + 0.818422*t + 0.1*(a^2) + 1.82743*a*t +
        0.180933*b*t - 0.822644*t^2 + 1.0201*a^2*t + 0.202*a*b*t + 0.01*(b^2)*t -
        0.74654*a*(t^2) + 0.82278*b*(t^2) + 0.522429*(t^3)
]
# The initial valuations DO contain floating point numbers, so
# approximate equality (≈) must be used against the expected results.
@assert isapprox(p_lie, actual, atol=0.00001)
@assert isapprox(p_pic, actual, atol=0.00001)



#
# Test the Lie derivative routine AND picard operator routine.
#
# Examine the outputs given an ODE that depends on time.

k = 3
set_variables("y t", order=2*k)
vars = get_variables(k)
y, t = vars


# The vector field contains goniometric functions AND depends on t!
#       y' = -y - sin(t) + cos(t)
# Now construct the taylorized vector field.
fpoly = [ -y - sin(t) + cos(t) ]

println("fpoly: ")
println(fpoly)

# The initial variable valuation; the bootstrapping values x0.
val0 = [ y ]

# Compute the Taylor polynomial approx. of the true flow.
p_lie = tay_poly(fpoly, k, vars, val0)
p_pic = picard(fpoly, k, vars, val0)

actual = [
    y + t - y*t - t^2 + (1/2)*y*t^2 + (1/6)*t^3,
]
# The initial valuations do NOT contain floating point numbers, so
# exact equality (==) is even possible against the expected results.
@assert p_lie == actual
@assert p_pic == actual

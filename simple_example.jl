using TaylorSeries  # Uses IntervalArithmetic
using TaylorModels  # Version req'd: https://github.com/gaperez64/taylormodels.jl

# Dynamics
g(y,t) = -y - sin(t) + cos(t)

# Initial state bounds
# y(0) = [1, 1]
# t(0) = [0, 0]

# Taylor variables (from TaylorSeries library)
y, t = set_variables("y t", order=4)

# We want to have a polynomial approximation of the dynamics centered around
# the initial values. The library centers around 0 by default, so for y we
# feed the function y -> y + 1.
g_poly_approx = g(y + 1, t)
print(g_poly_approx)





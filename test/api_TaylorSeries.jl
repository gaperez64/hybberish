using TaylorSeries

#
# Setup.
#
# Make sure Big-O is printed, which represents a each TaylorN's order member.
displayBigO(true)

#
# Test `TaylorSeries/parameters.jl` functions. This relates to object
# `_params_TaylorN_` of type `ParamsTaylorN` which is TaylorSeries'
# internal storage for `TaylorN` information.
#

# Examine the structure of the TaylorN internal parameters.
const info = TaylorSeries._params_TaylorN_
# The default order for constructed TaylorN objects.
@assert(info.order == get_order())
# The number of variables in use.
@assert(info.num_vars == get_numvars())
# The ordered list of string labels of the variables.
@assert(info.variable_names == get_variable_names())
# The ordered list of Symbols of the variables.
@assert(info.variable_symbols == get_variable_symbols())
println("Default values of `TaylorSeries._params_TaylorN_`:")
show_params_TaylorN()
println()

# TaylorSeries.jl internally stores which variables are defined.
# We compare the univariate (x) and multivariate (x, y, ...) cases.
usednames = ["x", "y"]
order = 3
for i=1:2
    varnames = usednames[1:i]
    println("Set TaylorSeries variables: " * join(varnames, ", "))
    vars = set_variables(join(varnames, " "), order=order)
    # Even if only one variable is set, the variables are still
    # of the TaylorN type.
    @assert(all(var -> typeof(var) <: TaylorN, vars))
    println("var types: $(map(var -> string(typeof(var)), vars))")


    # TaylorSeries defines parameter-less functions allow
    # us to access the internally stored variable information.
    @assert(varnames == get_variable_names())
    @assert(length(varnames) == get_numvars())
    @assert(vars == get_variables())
    println("var names: $(varnames)")
    println("var count: $(length(varnames))")
    println("variables: $(vars)")
    println()
end



#
# Subsequent calls to `set_variables` do not modify existing TaylorN objects.
#

xa, = set_variables("x", order=10)
xb, = set_variables("x", order=xa.order*2)
@assert(2*xa.order == xb.order)
println("xa.order = $(xa.order), xb.order = $(xb.order)\n")



#
# Test `TaylorSeries.fixorder(a, b)` manually.
# The TaylorSeries arithmetic uses it.
#

x, y = set_variables("x y", order=3)

# If both TaylorN objects have the same order, they remain unchanged.
@assert(x.order == y.order)
xres, yres = TaylorSeries.fixorder(x, y)
@assert(x == xres && y == yres)
println("fixorder($x, $y)")
println("      = ($xres, $yres)")

# The order of the objects differs, lower both of them to be the same order.
# This may result in truncation.
# Define two polynomials:
#   a = x     of order 2
#   b = y^3   of order 3
xlow = get_variables(2)[1] # Lower x's order
ypow3 = y^3
# We define them to have different order.
@assert(xlow.order != ypow3.order)
xres, yres = TaylorSeries.fixorder(xlow, ypow3)
# Now the orders should be equalized.
@assert(xres.order == yres.order)
# Lowering the order of polynomial b to 2 truncates its
# singular term "y^3", so b = 0 because no terms remain.
@assert(xlow == xres && yres == 0.0)
println("fixorder($xlow, $ypow3)")
println("      = ($xres, $yres)\n")



#
# Test `TaylorN` construction.
#

vars = set_variables("x y", order=5)
x5, y5 = vars

# You can specify the type of the polynomial coefficients.
@assert(typeof(get_variables()[1]) == TaylorN{Float64}) # default type
@assert(typeof(get_variables(Float32)[1]) == TaylorN{Float32})

# You can retrieve variables of a different order than that
# stored internally in TaylorSeries, without modifying that
# stored order.
x3, y3 = get_variables(3)
@assert(x3.order == 3 && y3.order == 3)
@assert(x5.order == 5 && y5.order == 5)

# You cannot obtain variables of order > the truncation order.
# `get_order()` yields the maximal term order for all TaylorN polynomials.
excessive_order = get_order() + 2
try
    get_variables(excessive_order)
    @assert(false) # Should not reach here!
catch
    println("The internally stored max order is $(get_order()), "*
            "cannot construct a TaylorN of order $excessive_order\n")
end

# There are shortcut functions to construct a TaylorN for the i-th variable.
x, = set_variables("x", order=3) # The usual way, but has side-effects!
xs = TaylorN(1) # A shorthand for TaylorN of var x, without side-effects.
# These TaylorN objects should be identical.
@assert(x == xs && x.order == xs.order)
println("`set_variables(...)` and `TaylorN(::Int)` create identical variable objects.")

println()



#
# Test `TaylorN` same-order arithmetic..
#
# TODO: Check dynamics of combining different-order polynomials.

x, y = set_variables("x y", order=10)
x1, y1 = get_variables(1)

# Test arithmetic with constants.
# addition or subtraction
res = x + 4
@assert(res.order == x.order && res[0] == 4) # res[0] is the constant term.
println("$x + 4 = $res")
res = -x
@assert(res.order == x.order && res[1][1] == -1) # res[1][1] is the coeff of x^1
println("-$x = $res")
# multiplication
res = y * 4
@assert(res.order == y.order && res[1][2] == 4) # res[1][2] is the coeff of y^1
println("$x * 4 = $res")
# substitution
res = (2*x - 4y)([1, 12])
@assert(res == (2 - 48))
println("(2*x - 4y)([1, 12]) = $res")

# Test order changes
# Same order terms, everything is preserved, the order stays 10.
res = x + y
@assert(res.order == x.order)
println("$x + $y = $res")
# Different order terms, order is lowered to 1!
res = x + y1
@assert(res.order == y1.order)
println("$x + $y1 = $res")
# Same for substitution, y^5 is truncated because the order becomes 1!
res = (x + y)([x1, y1])
@assert(res.order == y1.order)
println("($x + $y)([$x1, $y1]) = $res")

# Test evaluations
# Evaluate using only number values.
res = x([1, 1])
@assert(res == 1)
println("$x([1, 1]) = 1")
# Evaluate using a combination of number and TaylorN values.
res = x([1, x])
@assert(res == 1)
println("$x([1, $x]) = 1")
# Evaluate using only TaylorN values.
scalarN = TaylorN(7, get_order())  # Get the scalar "7" as a TaylorN
res = (x + y)([x, scalarN])  # This is fine, scalarN has enough dimensions
@assert(res.order == x.order && res[0] == 7 && res[1][1] == 1) # res[1][1] is the coeff of x^1
println("($x + $y)([$x, $scalarN]) = $res")
try
    x([x, 1])    # This is not fine, "1" has too few dimensions?
    @assert(false) # Should never reach here!
catch
end

# Evaluate the polynomial "x" at x=0.
# The two calling conventions should be equivalent;
# for TaylorN p: evaluate(p) == p()
@assert(evaluate(x) == x() == 0)
println("Two equivalent evaluation calling conventions: evaluate(x) == x() == 0\n")

# Shift the expansion point of a polynomial.
p = x + 2*y
pc = deepcopy(p)
update!(pc, [1, -1])
#=
   (x + 2y)([1, -1])
 = (x + 1) + 2*(y - 1)
 = (1 - 2) + x + 2y
 = -1 + x + 2y
  The updated polynomial only differs in the constant term!
=#
@assert(p[0] == 0 && pc[0] == -1 && p[1:end] == pc[1:end])
println("pre-shift:  $p")
println("post-shift: $pc\n")

# Test a complicated evaluation with different input orders.
x5, y5, = set_variables("x y z w", order=5)
z2 = get_variables(2)[3]
w3 = get_variables(3)[4]
cte0 = TaylorN(0, get_order()) # scalar 0 as a TaylorN
res = nothing
try
    # An error occurs in the `power_by_squaring!` function of TaylorSeries.
    # It expects all of its inputs to be of the same order, but this example
    # uses different orders!
    global res = (2*x5 - x5*y5 + y5^3)([1+z2, w3^2, cte0, cte0])
    @assert(false) # Should never reach here!
catch
    println("Evaluation with different input orders failed!")
end
z2 = get_variables()[3]
w3 = get_variables()[4]
res = (2*x5 - x5*y5 + y5^3)([1+z2, w3^2, cte0, cte0])
# TaylorN arithmetic propagates the lowest order at each operator.
@assert(res.order == minimum(get_order.([x5, y5, z2, w3])))
println("Evaluation with same input orders succeeded")
println("$(2*x5 - x5*y5 + y5^3)($([1+z2, w3^2, cte0, cte0])) = $res")


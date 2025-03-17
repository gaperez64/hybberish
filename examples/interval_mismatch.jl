using TaylorSeries
using TaylorModels
#using IntervalArithmetic
#using Polynomials


function example_2210_function(x)
    #return x - x^2
    return x * (interval(1, 1) - x)
end

x, = set_variables("x", order=4)
example_2210_polynom = example_2210_function(x)
#Polynomial([interval(0), interval(1,1), -interval(1,1)])

println("Consider the polynomial x-x^2 from example 2.2.10. from the Xin Chen thesis.")
println("Polynomial expr:", example_2210_polynom)

println("Eval func: ", example_2210_function(interval(-1, 1)))
println("Eval polynom: ", example_2210_polynom([interval(-1, 1)]))




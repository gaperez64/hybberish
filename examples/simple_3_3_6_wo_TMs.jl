using TaylorSeries
using TaylorModels
#using IntervalArithmetic
#using Polynomials

x, y, t = set_variables("x y t", order=4)
doms = [interval(-1,1), interval(-0.5,0.5), interval(0,0.02)]

p1 = 1 + y - x^2 * t - x * t^2 - (1/3)*t^3
p1e = -x^2 * t - x * t^2 - (1/3)*t^3
intp1e = p1e(doms)
println((intp1e + interval(-0.1,0.1)) * interval(0,0.02))

p2e = -2*t^2*y - 2*x*t*y
intp2e = p2e(doms)
println((intp2e + interval(-0.1311,0.1311)) * interval(0,0.02))


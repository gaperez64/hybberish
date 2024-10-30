using TaylorSeries  # IntervalArithmetic
using TaylorModels

# Assumes the last variable is t
function tay_poly(f, k)
    vars = get_variables()
    t = vars[end]
    # Let's pad f with a 1 at the end for t
    fp1 = copy(f)
    push!(fp1, 1)
    # Also, prepare a valuation vector with t=0
    # FIXME: We're cheating to obtain a zero with same order as the other
    # variables
    val0 = copy(vars)
    val0[end] = t - t
    # Prepare a first lie derivative and the result
    g = copy(vars)
    deleteat!(g, length(vars))
    # Start a vector function for the result
    res = copy(g)
    
    for i = 1:k
        println("g = $g")
        g = TaylorSeries.jacobian(g, vars) * fp1
        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
        res += term
    end
    
    return res
end


function tay_model(p, domain, k)
end

# Example 3.3.6
domain = IntervalBox([-1..1,-0.5..0.5,0..0.02])
zd = zero(domain)  # creates a hyperrectangle of 0's with same dimension
x, y, t = set_variables("x y t", order=3)
# focusing on the second dimension of the RHS of dynamics
f(a) = -a^2
# focusing on the second dimension of polynomial approx
p(a,b,c) = a + c + b * c
tm3 = TaylorModelN(p(x,y,t), -0.1..0.1, zd, domain)
# FIXME: this is needed as a hack to allow TM multiplication
# we update the order two twice the original one (or more)
x, y, t = set_variables("x y t", order=8)
ftm3 = f(tm3)
# I happen to know the total degree of the result is 3
# so we just need to get rid of that one, i.e. it is p_e
(polynomial(ftm3)[3](-1..1,-0.5..0.5,0..0.02) + remainder(ftm3)) * (0..0.02)

vars = set_variables("x y t", order=4)
f = [1 + vars[2], -vars[1]^2]
p = tay_poly(f, 4)
println("poly from new TM = $p")

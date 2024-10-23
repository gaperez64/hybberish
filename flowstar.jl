using TaylorSeries

# Assumes the last variable is t
function tay_poly(f, k)
    vars = get_variables()
    # Let's pad f with a 1 at the end for t
    fp1 = copy(f)
    push!(fp1, 1)
    # Also, prepare a valuation vector with t=0
    # FIXME: We're cheating to obtain a zero with same order as the other
    # variables
    val0 = copy(vars)
    val0[end] = vars[1] - vars[1]
    # Prepare a first lie derivative and the result
    g = copy(vars)
    deleteat!(g, length(vars))
    # Start a vector function for the result
    res = copy(g)
    
    for i = 1:k
        println("g = $g, vars = $vars")
        g = TaylorSeries.jacobian(g, vars) * fp1
        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
        println("term = $term")
        res += term
    end
    
    return res
end

vars = set_variables("x y t", order=5)
f = [1 + vars[2], -vars[1]^2]
p = tay_poly(f, 5)
println("poly from new TM = $p")

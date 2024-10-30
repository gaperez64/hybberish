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
        g = TaylorSeries.jacobian(g, vars) * fp1
        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
        res += term
    end
    
    return res
end


function tay_model_error(f, p, domain, k)
    zd = zero(domain)
    error = -0.1..0.1
    # to account for composition, we'll square the current order
    old_order = get_order()
    sqd_order = old_order * get_order(p)
    x, = set_variables("x", order=sqd_order)

    # TODO: start loop
    # construct the TM now
    tm = map((h) -> TaylorModelN(h, error, zd, domain), p)
    ftm = map((h) -> f(h), tm)
    # TODO: truncate and compute error
    # TODO: end loop after contraction check
    # TODO: add loop to continue contracting for a while

    # we now restore the order
    x, = set_variables("x", order=old_order)
    return error
end

#TODO: Clean below and call the above
#
# Example 3.3.6
domain = IntervalBox([-1..1,-0.5..0.5,0..0.02])
# I happen to know the total degree of the result is 3
# so we just need to get rid of that one, i.e. it is p_e
(polynomial(ftm3)[3](-1..1,-0.5..0.5,0..0.02) + remainder(ftm3)) * (0..0.02)

vars = set_variables("x y t", order=4)
f = [1 + vars[2], -vars[1]^2]
p = tay_poly(f, 4)
println("poly from new TM = $p")

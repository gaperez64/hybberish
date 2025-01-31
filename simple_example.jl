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
        println("Lie derivative:")
        println(map((h) -> evaluate(h, vars), g))
        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
        res += term
    end
    return res
end


"""Compute the remainder of the TM extension of the picard operator:
    Pf((p, I)) = J.
"""
function picard_tm_extension(vector_field_tms::Vector{T}, tmv, domain, k) where T <: TaylorModelN
    f::Vector{TaylorN} = polynomial.(vector_field_tms)
    errors::Vector{Interval} = remainder.(vector_field_tms)
    # FIXME: Is a copy needed? The substitution f(g(x, t), t)
    # seems to have side effects?
    tmv = copy(tmv)
    # Perform the substitution operation of the Picard operator:
    #   f(g(x, s), s).
    # FIXME: Hack to allow for higher degree terms in the intermediate computation
    y, t = set_variables("y t", order=k*2)
    ftm = map((fj) -> evaluate(fj, tmv), f)
    y, t = set_variables("y t", order=k)
    # Add the error of aproximating the true vector field.
    add_to_rem(tm, err) = TaylorModelN(polynomial(tm), remainder(tm) + err, tm.x0, domain)
    ftm = add_to_rem.(ftm, errors)
    # Compute the integral's remainder part:
    #   (int_enclosure(pe) + I) * [0, t].
    # FIXME: Once again we assume the t interval is the last one
    tdom = domain.v[end]
    # Prepare two functions to define new error intervals pointwise
    int_enclosure(p_error_terms, dom) = sum(perr(dom) for perr in p_error_terms)
    pe(tm::TaylorModelN) = polynomial(tm)[k:end]
    return map((tmj) ->
	       (int_enclosure(pe(tmj), domain) + remainder(tmj)) * tdom, ftm)
end


# Dynamics
f_dot(y, t) = -y - sin(t) + cos(t)

# Integration task specification
# a. Take delta_t = 1
# b. Construct a flowpipe consisting of 4 Taylor models, including the initial
#    one
# c. Work with order 4
ord = 4
domy = -2..2

# Taylor variables (from TaylorSeries library)
y, t = set_variables("y t", order=ord)

# Initial state variable bounds and domain
# y(0) = [1, 1]
# t(0) = [0, 0]
vals = IntervalBox(1..1, 0..0)
doms = IntervalBox(domy, 0..0.1)

# Iteration 1
#
# Step 0: Taylorize the dynamics
# We want to have a polynomial approximation of the dynamics centered around
# the midpoint of the current values.
ytm = TaylorModelN(y, interval(0), vals, doms)
ttm = TaylorModelN(t, interval(0), vals, doms)
# FIXME: Hack to allow for higher degree terms in the intermediate computation
y, t = set_variables("y t", order=ord*2)
ftm = f_dot(ytm, ttm)
# FIXME: We go back to the lower degree afterwards
y, t = set_variables("y t", order=ord)
println("taylorized vector field/dynamics:")
println(ftm)

# Step 1: Obtain the polynomial part of the Taylor model
p = tay_poly([polynomial(ftm)], ord)
println("polynomial part of TM:")
println(p)

# Step 2: Obtain the remainder/error interval of the TM
# Start Picard iteration, we need the candidate/guessed TM
candidate_tm = TaylorModelN(p[1], -0.1..0.1, vals, doms)
# Then we take the TM extension of the approx'd vector field composed with the
# candidate TM
rems = picard_tm_extension([ftm], [candidate_tm, ttm], doms, ord)
println("error interval part of TM:")
println(rems)
if all(issubset.(rems, [remainder(candidate_tm)]))
    println("Contractive!")
end

# Step 3: Get the new local values (and interval box) and update domain for next step
# i.e. just change the domain of the time variable in doms
# TODO

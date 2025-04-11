"""This file contains examples for our TM integration algorithm helper functions.
"""

include("../src/tm_integration.jl")



#
# Test TM integration TaylorModelN identity map auxiliary function.
#

k = 5
set_variables("x y", order=2*k)
vars = get_variables(k)
x, y = vars

doms = IntervalBox(-105.5..105.5, -402..402)

# Note that Taylor model arithmetic may throw an `iscontained` assertion error.
# This indicates that the evaluation input values are not contained inside the
# domain of the evaluated TaylorModelN.
# You MUST choose values that fall inside the domain of the idmap TaylorModelN!
box = IntervalBox(-3.3..1.3, 7..12)
num = [21.486, -88]

# The identity map Taylor model should acts as the
# identity function `id(x) = x` with respect to
# Taylor model evaluation using any value vector.
idvec = id(vars, doms)

# The identity map should work for a Number vector input.
res = [idj(num) for idj in idvec]
@assert res[1] == num[1] && res[2] == num[2]

# The identity map should work for a IntervalBox input.
res = [idj(box) for idj in idvec]
@assert res[1] == box[1] && res[2] == box[2]

# The identity map should work for a TaylorModelN input.
res = [idj(idvec) for idj in idvec]
@assert res[1] == idvec[1] && res[2] == idvec[2]



#
# Test TM integration TaylorModelN matrix scaling.
#

k = 5
set_variables("x y", order=2*k)
vars = get_variables(k)
x, y = vars

doms = IntervalBox((-0.0657183)..(0.0657890), (-0.0488922)..(0.0491184))
tmx = TaylorModelN(  0.904667 + 0.982762*x - 0.184876*y, -0.5..0.5, doms)
tmy = TaylorModelN( -0.909333 + 0.184876*x + 0.982762*y, -2..4, doms)
tmv = [ tmx, tmy ]

# Scale the Taylot models to bound their ranges inside the box [-1, 1]^2.
scaling = scale(tmv)
computed = scaling * tmv
computed = [tm() for tm in computed]

# A scaling matrix is a diagonal matrix by definition.
@assert isdiag(scaling) "The scaling matrix is not a diagonal matrix."
# The range of each of the scaled Taylor models should be contained in [-1, 1].
@assert all(c -> issubset(c, -1..1), computed) "Must be subset of [-1, 1]."



#
# Test TM integration TaylorModelN linear map construction.
#

k = 5
set_variables("x y", order=2*k)
vars = get_variables(k)
x, y = vars

# A matrix that represents a linear map.
lin_map::Matrix{Float64} = [
    4.0 -8.0;
    0.0  1.2
]

dom = IntervalBox(-2..2, 4..10)

# Construct a vector of TaylorModelN that are the equivalent representation
# of the linear map Matrix.
tm_lin_map = linear_map(lin_map, dom, vars)
polys_lin_map = polynomial.(tm_lin_map)

display(polys_lin_map)

# Given linear map (p, I).
# All remainders should be zero: I = [0, 0].
@assert all(rem -> rem == 0..0, remainder.(tm_lin_map))
# The linear map may NOT have order 0 terms.
@assert all(p -> p == 0, constant_polynomial(polys_lin_map))
# The linear map may NOT have order > 1 terms.
@assert all(p -> p == 0, nonlinear_polynomial(polys_lin_map))
# The linear map must have the expected linear terms.
@assert polys_lin_map[1] == (4.0*x - 8.0*y)
@assert polys_lin_map[2] == (0.0*x + 1.2*y)



#
# The QR preconditioning example from the M. Neher 2006 Taylor models paper.
#

k=5
set_variables("x y t", order=2*k)
vars = get_variables(k)
x, y, t = vars

# The time domain should not matter for preconditioning.
# The inputs to preconditioning should have fixed the time variable t=t_(n+1),
# to that time variable t does not appear in the polynomial parts of TMs.
bogus_time_domain = (-Inf)..(Inf)

# Center the domain on zero.
# We need to counterbalance this shift when constructing the left initial TMs.
initial_dom = IntervalBox(0.95..1.05, -1.05..(-0.95), bogus_time_domain)
initial_mid = mid(initial_dom)
initial_dom = IntervalBox((initial_dom.v - initial_mid)...)
@assert length(vars) == length(initial_dom) "Length mismatch."
shifted_vars = [v+i for (v,i) in zip(vars, initial_mid)]

# The right TMs are an identity map in the space (ODE) variables.
Ur0 = id(vars[1:end-1], initial_dom)
# The left TMs are the Taylor model representation of the interval initial set.
Ul0 = id(shifted_vars[1:end-1], initial_dom)

println("### Initial TM vectors.")
println("Ul ="); display(Ul0); println()
println("Ur ="); display(Ur0); println()
println("domain(Ul) ="); display(domain.(Ul0)); println()
println("domain(Ur) ="); display(domain.(Ur0)); println()



rng = IntervalBox([tm() for tm in Ur0]...)
# Explicitly enforce equal domains so we can pick the domain of a TM.
@assert allequal(domain.(Ur0)) "All domains must be equal."
# Verify that the composition of Ul0 and Ur0 is validly preconditioned.
# The composition is called preconditioned iff, `Rng(Ur0) ⊆ domain(Ul0)`.
# FIXME: Cite definition 5.2 in M. NEHER.
dom = IntervalBox(domain(Ul0[1]).v[1:2]...)
@assert issubset(rng, dom) "$rng SHOULD SUBSETEQ $(dom)"



# The integrated left taylor model, the result of naive TM integration.
# These Taylor models should be the ones from the example above Algorithm 6.1
# in the paper by Neher, M. et al.
dom = IntervalBox(-0.05..0.05, -0.05..0.05, bogus_time_domain)
Ul1 = [
    TaylorModelN(
         0.904667 + 1.01*x + 0.10*y,
         (-5.09307e-5)..(7.86167e-5),
         dom),
    TaylorModelN(
        -0.909333 + 0.19*x + 1.01*y + 0.1*x^2,
        (-1.75707e-4)..(1.60933e-4),
        dom),
]
_, Ul1, Ur1 = precondition(Ul1, Ur0, vars)


println("#### Preconditioning output ITERATION 1")
println("domain(Ul1) ="); display(domain.(Ul1)); println()
println("Ul1 ="); display(Ul1); println()
println("range(Ul1) ="); display([tm() for tm in Ul1]); println()
println("domain(Ur1) ="); display(domain.(Ur1)); println()
println("Ur1 ="); display(Ur1); println()
println("range(Ur1) ="); display([tm() for tm in Ur1]); println()


# The integrated left taylor model, the result of naive TM integration.
# These Taylor models should be the ones from the example above Algorithm 6.1
# in the paper by Neher, M. et al.
dom = IntervalBox(-1..1, -1..1, bogus_time_domain)
Ul2 = [
    TaylorModelN(
         0.817551 + 0.0664561*x - 0.00433580*y,
         (-5.72276e-5)..(9.15947e-5),
         dom),
    TaylorModelN(
        -0.835195 + 0.0233831*x + 0.0471479*y +0.000418026*x^2 -
            0.000117424*x*y + 0.00000824612*y^2,
        (-1.80914e-4)..(1.80850e-4),
        dom),
]
_, Ul2, Ur2 = precondition(Ul2, Ur1, vars)

println("#### Preconditioning output ITERATION 2")
println("domain(Ul2) ="); display(domain.(Ul2)); println()
println("Ul2 ="); display(Ul2); println()
println("range(Ul2) ="); display([tm() for tm in Ul2]); println()
println("domain(Ur2) ="); display(domain.(Ur2)); println()
println("Ur2 ="); display(Ur2); println()
println("range(Ur2) ="); display([tm() for tm in Ur2]); println()

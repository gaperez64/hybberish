"""This file contains examples for our TM integration algorithm helper functions.
"""

include("../src/tm_integration.jl")


"""Normalize the polynomials via affine transformation so that the domains become [-1, 1]^n."""
normalize_taylor(tmv::Vector{TaylorModelN{N,T,S}}) where {N,T,S} = [
    TaylorModelN(
        TaylorSeries.normalize_taylor(polynomial(tm), domain(tm)),
        remainder(tm),
        IntervalBox(fill(-1..1, length(domain(tm)))...)
    )
    for tm in tmv
]

"""Construct the TM representation of a Matrix linear map.

    A linear map Matrix is equivalently a linear Taylor model (p, I) with
    the zero constant part and zero remainder. In other words, p is only
    allowed to contain terms of exactly order 1, and I = [0, 0].
"""
# Instead of an explicit order value, make use of the variable TaylorN
# objects to implicitly construct a TaylorN of the correct order.
# The vars MUST be the variable objects currently in use.
function linear_map(linear_coeffs::Matrix{T}, dom::IntervalBox{N,S}, vars::Vector{TaylorN{T}}) where {N,T,S}
    @assert length(vars) == length(dom) == size(linear_coeffs)[1]
    return [
        TaylorModelN(
            dot(coeffs, vars), 0..0, dom)
        for coeffs::Vector{T} in eachrow(linear_coeffs)
    ]
end



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
# TODO: Implement algo 6.1 part 2.
#

"""Perform step 2 of Algorithm 6.1 (QR Preconditioned Taylor model method).

    This routine follows section 3.3.2 of Xin Chen's thesis, and makes use of
    the references cited there by Xin Chen.
        1) [Loh92] R. J. Lohner. Computation of guaranteed enclosures for the
                   solutions of ordinary initial and boundary value problems.
                   In J. R. Cash et al., editor, Computational ordinary
                   differential equations, pages 425-435. Clarendon Press,
                   1992.
        2) [MB05]  K. Makino and M. Berz. Suppression of the wrapping effect by
                   taylor modelbased verified integrators: Long-term
                   stabilization by preconditioning. International Journal of
                   Differential Equations and Applications, 10(4):353-384,
                   2005.


    Also see "6. Preconditioned Quadratic Example." and
    "Algorithm 6.1 (QR preconditioned Taylor model method)"in the paper:
        Neher, M. & Jackson, Kenneth & Nedialkov, Nedialko. (2007).
        On Taylor model based integration of ODEs. SIAM J. Numerical Analysis.
        45. 236-262. 10.1137/050638448.


    Additionally, we take inspiration from the preconditioning in Flowstar.
    See the TaylorModelVec class constructor, which implements preconditioning.
    https://github.com/chenxin415/flowstar/blob/b85a3211748cb77b736fe4ad42ee02d8d2b81148/flowstar-toolbox/TaylorModel.h#L2838

        TaylorModelVec<DATA_TYPE>::TaylorModelVec(const std::vector<Interval> & box, std::vector<Interval> & domain)


    There algorithm recomputes the left and right Taylor model vectors.
    It must ensure that the inclusion property holds; the true flow must
    be contained in the recomputed left Taylor model vector.
    There are two possible solutions to ensure the inclusion property.
        1) Bound the right Taylor models.
           We can recompute the domain of the space (ODE) variables so
           that the inclusion property holds.
        2) Modify the left and right Taylor models via a scaling matrix
           so that the inclusion property holds.


    inclusion of the true flow in the left Taylor models.
    If true, bound the domain of the right model to recompute
    the domain of the left mode.

    @param[in] tmv_left The left Taylor models that were integrated in step 1.
    @param[in] tmv_right The right Taylor models that were untouched in step 1.
    @param[in] vars The space (ODE) variable objects.
    @param[in] bound_range A boolean switch for choosing how to satisfy the
                           inclusion property. If true, then recompute the
                           variable domains. If false, then use a scaling matrix.
    @return The pair of preconditioned left and right Taylor model vectors.
"""
function precondition(
        tmv_left::Vector{TaylorModelN{N,T,S}},
        tmv_right::Vector{TaylorModelN{N,T,S}},
        vars::Vector{TaylorN{T}},
        bound_range::Bool) where {N,T,S}

    println("## Inputs")
    println("Uljp1 ="); display(tmv_left); println()
    println("domain(Uljp1) ="); display(domain.(tmv_left)); println()
    println("Urj ="); display(tmv_right); println()
    println("domain(Urj) ="); display(domain.(tmv_right)); println("\n\n")

    # We rely on unmodified variables to construct linear maps.
    @assert vars == get_variables() "The variables must be unmodified."
    @assert length(tmv_left) == length(tmv_right) "Vector length mismatch."

    # Require all domains to be equal within a TaylorModelN vector.
    # FIXME: `allequal` requires Julia 1.8+. Add an explicit Julia version requirement?
    @assert allequal(domain.(tmv_left)) "All left TM domains must be equal."
    @assert allequal(domain.(tmv_right)) "All right TM domains must be equal."

    # TODO: At what point should the normalization of domain occur?
    unit_box = IntervalBox(fill(-1..1, length(vars))...)
    @assert(all([domain(tm) != unit_box for tm in tmv_left]),
        "The left Taylor models must NOT have normalized domains as a convention.")
    @assert(all([domain(tm) == unit_box for tm in tmv_right]),
        "The right Taylor models MUST have normalized domain by definition.")

    tmv_left = normalize_taylor(tmv_left)

    println("## Normalized Inputs")
    println("LEFT ="); display(tmv_left)
    display(domain(tmv_left[1])); println()
    println("RIGHT ="); display(tmv_right)
    display(domain(tmv_right[1])); println("\n\n")

    # Note, Julia's broadcasting for vectors will be heavily used
    # explicitly (broadcasting operator) and implicitly (function that
    # operate on vectors) to avoid excessive for loops or list comprehension.
    #
    # We asserted that the domains within a TaylorModelN vec were the same.
    # So one domain object per vector suffices.
    dom_left::IntervalBox{N,S} = domain(tmv_left[1])


    # (I)
    # Compute the QR factorization of the linear part of U_{l,j+1}
    # Each matrix row corresponds to one of the polynomials' linear terms.
    # FIXME: `stack` requires Julia 1.9+. Add an explicit Julia version requirement?
    linear_coeffs::Matrix{T} = stack([p[1].coeffs for p in polynomial.(tmv_left)], dims=1)
    Q::Matrix{T} = Matrix(qr(linear_coeffs).Q)

    # TODO: Delete this!
    # The Q matrix computed in the example has the same as the once computed above,
    # but the signs of some of the numbers are different.
    Q = [0.982762 -0.184876; 0.184876 0.982762]

    println("## (I)")
    println("linear coeffs ="); display(linear_coeffs); println()
    println("Q ="); display(Q); println("\n\n")

    # (II)
    #
    # Shift (add) all non-constant terms AND remainders of U_{l,j+1} to U_{r,j}.
    tmv_shift  = tmv_left - constant_polynomial(tmv_left)
    tmv_right += tmv_shift
    tmv_left  -= tmv_shift
    # The remainders were shifted from left to right, so set left's to zero.
    # In interval arithmetic generally `[a, b] - [a, b] != [0, 0]`.
    tmv_left = [ TaylorModelN(tm, 0..0) for tm in tmv_left ]

    println("## (II.1)")
    println("LEFT  ="); display(tmv_left); println()
    println("RIGHT ="); display(tmv_right); println("\n\n")

    # Make Q the linear part of U_{l, j+1}.
    tmv_left += linear_map(Q, dom_left, vars)

    println("## (II.2)")
    println("LEFT ="); display(tmv_left); println("\n\n")

    # Apply Q^-1 on U_{r,j}.
    QT::Matrix{T} = transpose(Q)
    # FIXME: Use an arbitrary domain. See a fixme below for more info.
    #        The linear map's domain is NOT involved in remainder computation;
    #        This linear map is just a vehicle for TM evaluation.
    arbitrary_domain = IntervalBox([tm() for tm in tmv_right]...)
    QT_tmv = linear_map(QT, arbitrary_domain, vars)
    tmv_right = [ tm(tmv_right) for tm in QT_tmv ]

    println("## (II.3)")
    println("Q^T ="); display(QT); println()
    println("RIGHT ="); display(tmv_right); println("\n\n")

    if bound_range
        # (III)
        #
        # Bound the range of the new U_{r,j}.
        # In other words, modify the domain of the independent variables.
        dom_left = IntervalBox([tm() for tm in tmv_right]...)
        tmv_left = [ TaylorModelN(tm, dom_left) for tm in tmv_left ]

    else
        # (IV)
        # Apply a scaling matrix S_(j+1) on U_{r,j} such that each component of the
        # range of U_{r,j+1} := S_(j+1)^-1 ◦ U_{r,j} is contained in [-1, 1] and
        # spans [-1, 1] approximately.
        Sinv = scale(tmv_right)
        # Convert the matrix to an equivalent form: a linear map Taylor model vector.
        #= 
           FIXME: We really don't care about the domain this linear map!
           ==> To satisfy the requirement of the composition "g(f(x))" that
           the range of "f(x)" must be contained in the domain of "g", we sadly
           have to provide some arbitrary domain value that satisfies the
           requirement.
           ==> Cheat by making the range of the TMs domain of the linear map.
           The "iscontained" check will exactly check that the range is contained
           in the domain, so this choice of domain will always work. =#
        arbitrary_domain = IntervalBox([tm() for tm in tmv_right]...)
        Sinv_tmv = linear_map(Sinv, arbitrary_domain, vars)
        # Compute U_{r,j+1} := S ◦ U_{r,j} so that each component of the range of
        # the composition is contained in [-1, 1].
        # FIXME: But, does the range also span approximately [-1, 1]?
        tmv_right = [tm(tmv_right) for tm in Sinv_tmv]
        
        # (V)
        S_ = inv(Sinv)
        # Convert the matrix to an equivalent form: a linear map Taylor model vector.
        #=
          FIXME: Will this linear map domain ever trigger an "iscontained" assertion?
          The domain of this linear map matters. It will be composed into
          another Taylor model, so the map will be involved in Taylor model
          arithmetic and as such influences the computed remainders. =#
        S_tmv = linear_map(S_, unit_box, vars)
        # Set U_{l, j+1} := U_{l, j+1} ◦ S_(j+1)
        tmv_left = [tm(S_tmv) for tm in tmv_left]
    end

    return tmv_left, tmv_right
end



k = 10
set_variables("x y", order=2*k)
vars = get_variables(k)

# TODO: Construct U_{l,0} and U_{r,0}
# TODO Use  0.051  as bound instead of 0.O5 to trick the precondition
# assertion below. IntervalArithmetic is being annoying.
initial_dom = IntervalBox(-0.051..0.051, -0.051..0.051)
Ul0 = id(vars, initial_dom) # TODO: Is this how to construct Ul0????
Ul0 = scale(Ul0)*Ul0 # FIXME: Is it even allowed to scale the initial set???

# TODO: normalize (TaylorSeries.normalize_taylor) these TMs?
#       Because it should have domains [-1, 1]^n ???
Ur0 = id(vars, initial_dom)
Ur0 = normalize_taylor(Ur0)

rng = IntervalBox([tm() for tm in Ul0]...)
# Explicitly enforce equal domains so we can pick the domain of a TM.
@assert allequal(domain.(Ur0)) "All domains must be equal."
# Verify that the composition of Ul0 and Ur0 is validly preconditioned.
# The composition is called preconditioned iff, `Rng(Ul0) ⊆ domain(Ur0)`.
# FIXME: Cite definition 5.2 in M. NEHER.
@assert issubset(rng, domain(Ur0[1])) "$rng SUBEQ $(domain(Ur0[1]))"




# The integrated left taylor model, the result of naive TM integration.
# These Taylor models should be the ones from the example above Algorithm 6.1
# in the paper by Neher, M. et al.
dom = IntervalBox(-0.05..0.05, -0.05..0.05)
Uljp1 = [
    TaylorModelN(
         0.904667 + 1.01*x + 0.10*y,
         (-5.09307e-5)..(7.86167e-5),
         dom),
    TaylorModelN(
        -0.909333 + 0.19*x + 1.01*y + 0.1*x^2,
        (-1.75707e-4)..(1.60933e-4),
        dom),
]
# FIXME: Is this even a correct identity map?
#        Because the variables are just 0 constants!
Ur0 = id([0 * v for v in vars], IntervalBox(fill(-1..1, length(dom))...))

# TODO: For now ignore Ul0, just plug in the values from the example
# in the paper.
Uljp1, Urjp1 = precondition(Uljp1, Ur0, vars, false)

println("domain(Uljp1) ="); display(domain.(Uljp1)); println()
println("Uljp1 ="); display(Uljp1); println()
println("domain(Urjp1) ="); display(domain.(Urjp1)); println()
println("Urjp1 ="); display(Urjp1); println()

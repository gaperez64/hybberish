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



    @param[in] Uljp1 The left Taylor model that was integrated in step 1.
    @param[in] Urj The right Taylor model that was untouched in step 1.
"""
function precondition(
        Uljp1::Vector{TaylorModelN{N,T,S}},
        Urj::Vector{TaylorModelN{N,T,S}}) where {N,T,S}

    #= TODO: Note the following naming convention, used for brevity.

        Uljp1 stands for Ũ_{l, j+1} (note the tilde!).
        Urj stands for U_{r, j}
    =#

    # FIXME: Assume this is a vector of polynomial references.
    Uljp1_polys::Vector{TaylorN{T}} = [ polynomial(tm) for tm in Uljp1 ]
    Urj_polys::Vector{TaylorN{T}} = [ polynomial(tm) for tm in Urj ]


    # TODO: Consider using these methods.
    TaylorSeries.linear_polynomial
    TaylorSeries.nonlinear_polynomial
    TaylorSeries.constant_term

    # TODO: Rewrite the comments of the steps below as mathematical assignments
    # or mathematical operations.

    # (I)
    # Compute the QR factorization of the linear part of U_{l,j+1}
    # Each row of the matrix corresponds to one of the Taylor models.
    # FIXME: `stack` requires Julia 1.9+
    linear_coeffs::Matrix{T} = stack([p[1].coeffs for p in Uljp1_polys], dims=1)
    Q::Matrix{T}, _ = qr(linear_coeffs)

    # (II)
    # Shift all but the constant part of U_{l,j+1} to U_{r,j}.
    # TODO: Does shift mean 1) add to the existing poly OR 2) replace the existing poly???
    Urj_polys = [
        constant_term(pr) + linear_polynomial(pl) + nonlinear_polynomial(pl)
        for (pl, pr) in zip(Uljp1_polys, Urj_polys)
    ]
    # Make Q the linear part of U_{l, j+1}.
    # TODO: i.e. U_{l, j+1} := constant(U_{l, j+1}) + Q
    @assert length(Uljp1_polys) == size(Q)[1]
    Uljp1_polys = [
        TaylorN(constant_term(pl), get_order(pl)) +
        HomogeneousPolynomial(coeffs, 1)
        for (pl::TaylorN{T}, coeffs::Vector{T}) in zip(Uljp1_polys, eachrow(Q))
    ]
    # TODO: Apply Q^-1 on U_{r,j}.

    # (III)
    # TODO: Bound the range of the new U_{r,j}.

    # (IV)
    # Apply a scaling matrix S_(j+1) on U_{r,j} such that each component of the
    # range of U_{r,j+1} := S_(j+1)^-1 ◦ U_{r,j} is contained in [-1, 1] and
    # spans [-2, 1] approximately.

    # (V)
    # Set U_{l, j+1} := U_{l, j+1} ◦ S_(j+1)
end

x, y = set_variables("x y", order=10)
Uljp1 = [
    TaylorModelN(1.0 + 2.0*x - 3.0*y + x^2, 0..0, IntervalBox(-1..1, -2..2)),
    TaylorModelN(1.0 - 3.0*x + 2.0*y + y^2, 0..0, IntervalBox(-1..1, -2..2)),
]
Urj = [
    TaylorModelN(x, 0..0, IntervalBox(-1..1, -2..2)),
    TaylorModelN(y, 0..0, IntervalBox(-1..1, -2..2)),
]

precondition(Uljp1, Urj)

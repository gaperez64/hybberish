using LinearAlgebra

include("taylor_models/BasicTaylorModels.jl")


"""Generate the unit IntervalBox [-1, 1]^n where n = length(a)."""
unitbox(a) = IntervalBox(fill(-1..1, length(a)))

"""Normalize the polynomials via affine transformation so that the domains become [-1, 1]^n."""
normalize_taylor(tmv::Vector{TaylorModelN{N,T,S}}) where {N,T,S} = [
    TaylorModelN(
        TaylorSeries.normalize_taylor(polynomial(tm), domain(tm)),
        remainder(tm),
        unitbox(domain(tm))
    )
    for tm in tmv
]

"""Construct the TM representation of a Matrix linear map.

    A linear map Matrix is equivalently a linear Taylor model (p, I) with
    the zero constant part and zero remainder. In other words, p is only
    allowed to contain terms of exactly order 1, and I = [0, 0].

    Instead of an explicit order value, make use of the variable TaylorN
    objects to implicitly construct a TaylorN of the correct order.
    The vars MUST be the variable objects currently in use.
"""
function linear_map(linear_coeffs::Matrix{T}, dom::IntervalBox{N,S}, vars::Vector{TaylorN{T}}) where {N,T,S}
    @assert(length(vars) == size(linear_coeffs)[1],
        "Each linear expression must have one component for each var.")
    @assert(length(vars) <= length(dom),
        "Each var must define a domain component.")
    return [
        TaylorModelN(dot(coeffs, vars), 0..0, dom)
        for coeffs::Vector{T} in eachrow(linear_coeffs)
    ]
end

"""Create the identity map of an interval box.

    The identity map is a Taylor model vector (p, I) so that `p(box) + I = box`.
    In more general terms, (p, I) represents the identity function `id(x) = x`.
    Practically, the result is a vector of identity functions
        (p, I) = [
            (p_1, I_1),
            ...,
            (p_n, I_n)
        ]
    where `(p_j, I_j)([x_1, ..., x_n]) = x_j` is the j-th identity function.

    Note that we require the variable objects and domains to be passed.
    Due to constraints on TaylorModelN arithmetic, these must be known
    Explicitly at the time of construction of the identity map.

    @param[in] vars The ODE variable objects to construct the polynomial part.
    @param[in] doms The domains to assign the identity map Taylor models.
    @return The identity map.
"""
function id(vars::Vector{TaylorN{T}}, doms::IntervalBox{N,S})::Vector{TaylorModelN{N,T,S}} where {N,T,S}
    @assert(length(vars) <= length(doms),
        "Each variable must specify its domain for the identity map.")
    return [
        # Choose `(p, I) = (x, [0, 0])`  so that  `(p, I)(b) = b + [0, 0] = b`.
        # This is exactly the identity map `id(b) = b`.
        # This property should also hold for Taylor model composition.
        TaylorModelN(var, 0..0, doms)
        for var in vars
    ]
end

"""Compute a scaling matrix S that bounds `Rng(S*tmv)` within box [-1, 1]^n.

    Scaling matrix S is a diagonal matrix by definition. Suppose S is nxn=3x3.
            | s_1   0    0  |
        S = |  0   s_2   0  |
            |  0    0   s_3 |
    Then the scaling factor s_i is chosen so scaling the i-th component
    of the Taylor model vector by this factor bounds its range within [-1, 1],
    `Rng(s_i * tmv[i]) in [-1, 1]`.

    Note that `Rng(tmv[i]) = [0, 0]` would result in `s_i = 1/0 = Inf`.
    Since `[0, 0] in [-1, 1]` already holds, default to `s_i = 1` instead.
    All other degenerate interval ranges `Rng(tmv[i]) = [a, a]` where `a != 0`
    result in `s_i = 1/a != Inf` so that no ad hoc correction is needed.

    @param[in] tmv The Taylor models whose range to bound within [-1, 1].
    @return The scaling matrix S.
"""
function scale(tmv::Vector{TaylorModelN{N,T,S}})::Matrix{T} where {N,T,S}
    # The case where `Rng(tmv[i]) == [0, 0]` would result in `s_i = 1/0 = Inf`.
    # Since `[0, 0] in [-1, 1]` already, default to `s_i = 1` instead.
    nozero = v -> v==0 ? 1.0 : v
    # Suppose `Rng(tm_i) = [a, b]` then `s_i = 1 / max{ abs(a), abs(b) }`.
    # But, computing the scale value `s_i` requires evaluating a division 1/x.
    # This introduces a potential rounding error. To counteract this, to
    # ensure the range bounds do not become e.g. [-1.00001, 1.00001], make
    # `s_i` slightly larger to account for the rounding error.
    return diagm([ 1.0 / (nozero(mag(tm())) * 1.001) for tm in tmv ])
end

"""Perform step 2 of Algorithm 6.1 (QR Preconditioned Taylor model method).

    This routine draw from section 3.3.2 of Xin Chen's thesis, and makes use of
    the references cited there by Xin Chen.
        1) [MB05]  K. Makino and M. Berz. Suppression of the wrapping effect by
                   taylor modelbased verified integrators: Long-term
                   stabilization by preconditioning. International Journal of
                   Differential Equations and Applications, 10(4):353-384,
                   2005.

    Also see "6. Preconditioned Quadratic Example." and
    "Algorithm 6.1 (QR preconditioned Taylor model method)"in the paper:
        Neher, M. & Jackson, Kenneth & Nedialkov, Nedialko. (2007).
        On Taylor model based integration of ODEs. SIAM J. Numerical Analysis.
        45. 236-262. 10.1137/050638448.

    @param[in] tmv_left The left Taylor models that were integrated in step 1.
    @param[in] tmv_right The right Taylor models that were untouched in step 1.
    @param[in] vars The space (ODE) variable objects, plus the time variable object.
    @return The pair of preconditioned left and right Taylor model vectors.
"""
function precondition(
        tmv_left::Vector{TaylorModelN{N,T,S}},
        tmv_right::Vector{TaylorModelN{N,T,S}},
        vars::Vector{TaylorN{T}}) where {N,T,S}

    # We rely on unmodified variables to construct linear maps.
    @assert vars == get_variables() "The variables must be unmodified."
    @assert length(tmv_left) == length(tmv_right) "Vector length mismatch."

    # Require all domains to be equal within a TaylorModelN vector.
    # FIXME: `allequal` requires Julia 1.8+. Add an explicit Julia version requirement?
    @assert allequal(domain.(tmv_left)) "All left TM domains must be equal."
    @assert allequal(domain.(tmv_right)) "All right TM domains must be equal."

    # Assume the last variable represents time t.
    t = vars[end]
    vars_no_t = vars[1:end-1]
    # NOTE: We assume that neither the left nor right Taylor models contain
    # time variable t in the polynomial parts. So whatever the value of or
    # whatever the modifications made to the time domains, it does not affect
    # the intermediate or final values for preconditioning.
    timedoml = domain(tmv_left[1]).v[end]
    unitdom = unitbox(vars)

    # FIXME: Should normalization occur inside this function, or outside?
    tmv_left = normalize_taylor(tmv_left)
    # Construct a dummy TM because TM composition requires one TM per variable.
    ttm = TaylorModelN(t, 0..0, unitdom)

    # (I)
    # Compute the QR factorization of the linear part of U_{l,j+1}
    # Each matrix row corresponds to one of the polynomials' linear terms.
    # FIXME: `stack` requires Julia 1.9+. Add an explicit Julia version requirement?
    linear_coeffs::Matrix{T} = stack([p[1].coeffs for p in polynomial.(tmv_left)], dims=1)
    Q::Matrix{T} = Matrix(qr(linear_coeffs).Q)

    # FIXME: Ad hoc: modify the domain of tmv_right in anticipation that the
    # composition with the QT TMs will induce this domain.
    tmv_right = [ TaylorModelN(tm, unitdom) for tm in tmv_right ]

    # (II)
    #
    # Shift (add) all non-constant terms AND remainders of U_{l,j+1} to U_{r,j}.
    tmv_shift = tmv_left - constant_polynomial(tmv_left)
    tmv_right = [ tm(vcat(tmv_right, ttm)) for tm in tmv_shift ]
    tmv_left -= tmv_shift
    # The remainders were shifted from left to right, so set left's to zero.
    # In interval arithmetic generally `[a, b] - [a, b] != [0, 0]`.
    tmv_left = [ TaylorModelN(tm, 0..0) for tm in tmv_left ]

    # Make Q the linear part of U_{l, j+1}.
    tmv_left += linear_map(Q, unitdom, vars_no_t)

    # Apply Q^-1 on U_{r,j}.
    QT::Matrix{T} = transpose(Q)
    # FIXME: Use an arbitrary domain. See a fixme below for more info.
    #        The linear map's domain is NOT involved in remainder computation;
    #        This linear map is just a vehicle for TM evaluation.
    # FIXME: Don't use arbitrary domain! Use the domain of the left Taylor models
    # specifically, because Q originates from the left Taylor Models?
    arbitrary_domain = IntervalBox([tm() for tm in tmv_right]..., ttm())
    QT_tmv = linear_map(QT, arbitrary_domain, vars_no_t)
    tmv_right = [ tm(vcat(tmv_right, ttm)) for tm in QT_tmv ]

    # (III)
    #
    # Bound the range of the new U_{r,j}.
    # In other words, modify the domain of the independent variables.
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
    arbitrary_domain = IntervalBox([tm() for tm in tmv_right]..., ttm())
    Sinv_tmv = linear_map(Sinv, arbitrary_domain, vars_no_t)
    # Compute U_{r,j+1} := S ◦ U_{r,j} so that each component of the range of
    # the composition is contained in [-1, 1].
    # FIXME: But, does the range also span approximately [-1, 1]?
    tmv_right = [ tm(vcat(tmv_right, ttm)) for tm in Sinv_tmv ]

    # The scaled Taylor models should have a range subset of the unit box.
    # If this fails, then likely the computed scaling matrix does not properly
    # take rounding errors into account when computing the scale values.
    rng = IntervalBox([tm() for tm in tmv_right]...)
    dom = IntervalBox(unitdom.v[1:end-1]...)
    @assert issubset(rng, dom) "The scaled range is not subset of the unit box, got: $rng"

    # (V)
    S_ = inv(Sinv)
    # Convert the matrix to an equivalent form: a linear map Taylor model vector.
    #=
        FIXME: Will this linear map domain ever trigger an "iscontained" assertion?
        The domain of this linear map matters. It will be composed into
        another Taylor model, so the map will be involved in Taylor model
        arithmetic and as such influences the computed remainders. =#
    S_tmv = linear_map(S_, unitdom, vars_no_t)
    # Set U_{l, j+1} := U_{l, j+1} ◦ S_(j+1)
    tmv_left = [ tm(vcat(S_tmv, ttm)) for tm in tmv_left ]

    println("Range(tmv_right) ="); display(IntervalBox([tm() for tm in tmv_right]...)); println();

    # Retain the input time components.
    tmv_left  = [
        TaylorModelN(tm, IntervalBox(domain(tm).v[1:end-1]..., timedoml))
        for tm in tmv_left
    ]
    tmv_right = [
        TaylorModelN(tm, IntervalBox(domain(tm).v[1:end-1]..., timedoml))
        for tm in tmv_right
    ]

    # For practical purposes, the range of each component of the right
    # Taylor models should be contained in [-1, 1].
    # FIXME: Make the check approximate, we need to account for
    #        over-approximate arithmetic resulting in range bounds
    #        that are inside [-1, 1] with a small tolerance.
    @assert all(rng -> -1.01 < rng.lo && rng.hi < 1.01, [tm() for tm in tmv_right])

    return tmv_left, tmv_right
end

"""Generate the Taylor polynomial approximation of the true flow specified by
   given dynamics (ODEs) up to the given degree via Lie derivatives.

   We require TaylorN representations of the variables to be supplied by the
   caller. This allows the caller to decide their order; the variables' orders
   influence the orders of the output polynomials.

   @param[in]    f The Taylorized dynamics (ODEs).
   @param[in]    k The degree (order) of the resulting polynomials.
   @param[in] vars The TaylorN objects representing the variables.
        Assume the last element represents the time variable.
   @param[in] val0 The TaylorN objects representing the initial variable
        valuation. This vector contains NO time component.
   @return The polynomial approximations of the true flow.
"""
function tay_poly(
        f::Vector{TaylorN{N}},
        k::Integer,
        vars::Vector{TaylorN{N}},
        val0::Vector{TaylorN{N}}) where {N <: Number}
    vars_orders = get_order.(vars)

    # We require a vector of unmodified variables for derivation.
    @assert vars == get_variables()
    # The function `evaluate(::TaylorN, ::Vector{TaylorN})` does not play nice
    # with different orders for the substitution values. So, require all orders
    # to be the same.
    @assert all( vars_orders[1] .== vars_orders )
    # Taylor series arithmetic propagates the lowest order of its operands.
    # To generate order k polynomials, all variables must be at least order k.
    @assert all( k .<= vars_orders )
    # The variables contain an additional last component: the time variable.
    @assert length(f) == length(vars)-1
    @assert length(val0) == length(vars)-1

    t = vars[end]
    # Let's pad f with a 1 at the end for t
    fp1 = vcat(f, 1)
    # Prepare a first lie derivative and the result,
    # Note that it does not contain a time component.
    g = vars[1:end-1]
    # Start a vector function for the result.
    lie_derivative = copy(val0)
    # Also, prepare a valuation vector with t=0
    # Set time to a zero polynomial with the same order as other variables.
    val0 = vcat(val0, zero(t))
    for i = 1:k
        # The variables contain a time component,
        # so the jacobian contains a time column.
        g = TaylorSeries.jacobian(g, vars) * fp1
        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
        lie_derivative += term
    end
    println("Final Lie derivative:"); println(lie_derivative); println()
    return lie_derivative
end


"""Compute the remainder of the TM extension of the picard operator.

    The term "TM extension" refers to the fact that all mathematical operators
    in the Picard operator are replaced with their TM arithmetic counterparts.
    This implies that the vector field f and function g w.r.t. which the normal
    Picard operator is computed, must be TMs instead of functions when working
    with the TM extension of the Picard operator P.
    Suppose

        P_f(g) = x * int_0^t f(g(x, t), t) dt

    is the normal Picard operator. Then

        P_F(G)
        = x * int_0^t F(G) dt
        = x * int_0^t (p, I) \\circ (q, J) dt
        = x * int_0^t p(q, J) + I dt
        = x * int_0^t (r, K) + I dt
        = x * int_0^t (r, K + I) dt
        = (p, K + I)

    is the TM extension of the Picard operator, where F = (p, I) is a TM
    (Taylorization plus Lagrange remainder) of the vector field f and
    G = (q, J) is the TM representing the function g.

    @param[in] vector_field_tms The TM representation F of the vector field f.
    @param[in]     function_tms The TM representation G of the function g.
    @return The remainders of the TM result of the TM extension of the Picard
     operator; the interval vector K.
"""
function picard_tm_extension(
        p::Vector,
        J::IntervalBox{M,S},
        vars::Vector,
        doms::IntervalBox{N,S},
        Dli::Vector{TaylorModelN{N,T,S}},
        real_t::Float64,
        vector_field_constructor::Function) where {N,M,T,S}

    @assert length(J) == length(p)
    @assert length(J) == length(doms)-1

    # FIXME: Implement a general picard operator function so that it works
    #        both for constructing p given a TaylorN vector AND for doing
    #        the contractiveness check given a TaylorModelN vector.

    # Time-shift the time variable during Taylorization of time-dependent ODEs.
    # I.e. in the ODEs substitute variable t by Taylor model (t + ti) + [0, 0].
    candidate_ttm = TaylorModelN(vars[end] + real_t, interval(0), doms)
    # Compose the candidate Taylor models that we guess envelop the true flow.
    candidate_tmv = [ TaylorModelN(p_i, J_i, doms) for (p_i, J_i) in zip(p, J) ]

    #= Step (1), perform the composition operation of the
      TM extension of Picard operator, which  accounts for errors
      coming from the dynamics having been approximated by polynomials.
        (p, I) \circ (q, J) = p(q, J) + I = (r, K) + I = (r, K + I)
      The normal Picard operator expresses this as f(g(x, t), t).
    =#
    fgtmv = Vector{TaylorModelN{N,T,S}}(undef, length(p))
    vector_field_constructor(fgtmv, vcat(candidate_tmv, candidate_ttm))

    #= Step (2), apply the antiderivative formula.
        (Int(pe) + I) * \delta
      where `pe` represents the truncated terms, `Int(pe)` is its interval
      enclosure and `\delta` is the time step.
      Every TaylorModelN object `tm` specifies its own polynomial order `d`.
      This is functionally the truncation order of that Taylor model.
      Consequently, the integral operation would only make the
      HomogeneousPolynomial `pe = polynomial(tm)[end]` exceed the truncation order,
      since that polynomial contains all terms of order exactly `d`.
      If we truncate before we apply the integral, then `pe` corresponds exactly to
      the terms to truncate.
    =#
    """Truncated term interval enclosure for truncation before integration."""
    intpe(tm::TaylorModelN) = evaluate(polynomial(tm)[end], domain(tm))
    """Get the width of the time component of the TaylorModelN domain."""
    dt(tm::TaylorModelN) = diam(domain(tm).v[end])
    integral_rem = IntervalBox([
        ( intpe(tm) + remainder(tm) ) * dt(tm) for tm in fgtmv
    ])

    # The Picard operator depends on the current initial set.
    # Add the remainders of the current initial set to imitate
    # how the Picard operator would work.
    # FIXME: This should happen as part of a generalized Picard function.
    bootstrap_rem = IntervalBox(remainder.(Dli))

    return bootstrap_rem + integral_rem
end

"""Compute a safe remainder interval for the i-th flowpipe.

    The TM extension of the Picard operator, Pf, can be used to
    refine the remainder estimate of the i-th flowpipe.
    Suppose we compute the TM (p, J_{0}) before refinement.
    This function computes a safe remainder by first finding a
    contractive remainder
        (p, J_{m})  where m >= 0
    and then continually applying refinement
        Pf((pl, J_{j})) = (pl, J_{j+1})  for j = m, ..., n
    to finally obtain the contractive, refined TM
        (p, J_{n}).
    Note that m, n are usually not known beforehand, but instead
    become known when the contractiveness and refinement iteration
    converge or fail.

    @param[in] vector_field_constructor The vector field function.
    @param[in] p The true flow Taylor polynomial approximations.
    @param[in] J The initial remainder estimate.
    @param[in] vars The variable objects to use.
    @param[in] doms The previous interval inintial set, with updated time domain???
    @param[in] Dli The current initial set.
    @param[in] real_t The real time ti of the integration set, the time of the current initial set.
    @param[in] NR_CONTRACTIVENESS_TRIES The max number of contractiveness tries to attempt.
    @param[in] NR_REFINEMENTS The max number of refinements to perform.
    @param[in] REFINEMENT_EPS Quit refinement early if the improvement falls below this threshold.
    @param[in] SCALE The contractiveness widening scalar factor.
"""
function tay_model_error(
        vector_field_constructor::Function,
        p::Vector,
        J::Interval{S},
        vars::Vector,
        doms::IntervalBox{N,S},
        Dli::Vector{TaylorModelN{N,T,S}},
        real_t::Float64,
        NR_CONTRACTIVENESS_TRIES::Integer,
        NR_REFINEMENTS::Integer,
        REFINEMENT_EPS::Float64,
        SCALE::Float64) where {N, T, S}

    @assert NR_CONTRACTIVENESS_TRIES > 0
    @assert NR_REFINEMENTS >= 0
    @assert REFINEMENT_EPS > 0.0
    @assert SCALE > 1.0

    @assert length(p) == length(vars)-1

    # Setup.
    J0 = nothing
    J1 = nothing
    Jn = nothing

    # Start Picard iteration, we need the candidate/guessed TM
    for ctry in 1:NR_CONTRACTIVENESS_TRIES
        J0 = IntervalBox(fill(J, length(p))...)

        # The TM extension of the Picard operator reveals contractiveness
        # of the candidate safe Taylor model vector (p, J0).
        J1 = picard_tm_extension(p, J0, vars, doms, Dli, real_t, vector_field_constructor)

        # Test contractiveness.
        if all(issubset.(J1, J0))
            # Contractiveness test succeeded, pass along the contractive (safe)
            # remainder J1.
            break
        # Contractiveness failure condition reached: nr of retries exhausted.
        elseif ctry == NR_CONTRACTIVENESS_TRIES
            println("Could not find a contractive remainder.")
            @assert false "Could not find a contractive remainder after $ctry tries.\n"*
                          "J0 = $J0\nJ1 = $J1\n==> J1 is not a subset of J0."
        end

        J = J * SCALE
    end

    # Perform remainder refinement to tighten the bounds.
    Jn = J1
    for nr in 1:NR_REFINEMENTS
        print("Refinement no. $nr")
        Jprev = Jn

        # Refine (tighten) the safe remainder.
        Jn = picard_tm_extension(p, Jn, vars, doms, Dli, real_t, vector_field_constructor)

        @assert all(issubset.(Jn, Jprev)) "Refinement should only increase the bound tightness!"
        max_improvement::Float64 = maximum(diam.(Jprev) - diam.(Jn))

        println("  (max improvement=$max_improvement)")
        if max_improvement < REFINEMENT_EPS
            break
        end
    end

    return Jn
end

"""The QR preconditioned TM integration algorithm.

    Compute an overapproximation of the true flow of the
    system of ODEs represented by the given vector field,
    as a sequence of flowpipes over partial time horizons
    [0, δi] of the full time horizon [0, Δ].

    @param[in] vector_field_constructor A constructor for the TM vector
                                        representation of the vector field.
    @param[in]                  initial The names and domains of the ODE
                                        variables.
    @param[in]                        k The TM arithmetic/truncation order.
    @param[in]                        J The initial remainder estimate for
                                        finding a contractive remainder.
    @param[in]             time_horizon The bounded time horizon.
    @param[in]           TIME_STEP_SIZE The fixed time step size
    @param[in] NR_CONTRACTIVENESS_TRIES The number of times to attempt widening
                                        the remainder estimate in order to find
                                        a contractive remainder.
    @param[in]           NR_REFINEMENTS The number of additional remainder
                                        refinements to apply after finding the
                                        contractive remainder.
    @param[in]                    SCALE The scale factor by which to widen the
                                        remainder estimate when the
                                        contractiveness test fails.
    @param[in]           REFINEMENT_EPS The minimum improvement a remainder
                                        refinement step should affect. If the
                                        improvement falls below this threshold,
                                        then we declare the refinement to have
                                        converged and stop the refinment loop.
"""
function tm_integration_QR(
        vector_field_constructor::Function,
        initial::Vector{Tuple{String, Interval{S}}},
        k::Integer,
        J::Interval{S},
        TIME_HORIZON::Float64,
        TIME_STEP_SIZE::Float64,
        NR_CONTRACTIVENESS_TRIES::Integer,
        NR_REFINEMENTS::Integer;
        SCALE::Float64=2.0,
        REFINEMENT_EPS::Float64 = 0.001) where {S}

    @assert(TIME_HORIZON > 0)  # The time horizon must not be [0, 0].
    @assert(TIME_STEP_SIZE > 0)
    @assert(NR_CONTRACTIVENESS_TRIES >= 0)
    @assert(NR_REFINEMENTS >= 0)

    # Unpack the input parameters.
    numvars = length(initial)
    _names::NTuple{numvars, String},
    _init::NTuple{numvars, Interval{S}} = zip(initial...)
    names::String = join(_names, ' ')
    init::IntervalBox{numvars, S}  = IntervalBox(_init)

    # Double the truncation degree to obtain the TaylorSeries max order.
    # This accounts for order-related assertions applicable to Taylor
    # series arithmetic.
    # The added degrees effectively are a buffer to make TM arithmetic work.
    kbuffer = 2*k

    # Construct Taylor variables (from TaylorSeries library)
    # Fix the TaylorSeries internal, max order.
    set_variables(names, order=kbuffer)
    # Construct the variables with the actual truncation degree of choice.
    vars = get_variables(k)

    # By definition of algo 6.1 in M. Neher (2006), the symbolic (normalized)
    # space variables should have domain [-1, 1]^m.
    # So, normalize the variables to have domain [-1, 1]^m.
    # These variables take the form: c + x*A
    vars_normalized = [
        TaylorSeries.normalize_taylor(v, init)
        for v in vars[1:end-1]
    ]
    vars_normalized = vcat(vars_normalized, vars[end]) # Add the time variable.
    init = IntervalBox(unitbox(init.v[1:end-1])..., init.v[end])

    # The left TMs are the Taylor model representation of the interval initial set.
    Dl0 = id(vars_normalized[1:end-1], init)
    # The right TMs are an identity map in the space (ODE) variables.
    Dr0 = id(vars[1:end-1], init)

    rng = IntervalBox([tm() for tm in Dr0]...)
    # Explicitly enforce equal domains so we can pick the domain of a TM.
    @assert allequal(domain.(Dr0)) "All domains must be equal."
    # Verify that the composition of Dl0 and Dr0 is validly preconditioned.
    # The composition is called preconditioned iff, `Rng(Dr0) ⊆ domain(Dl0)`.
    # See definition 5.2 in M. NEHER (2006).
    dom = IntervalBox(domain(Dl0[1]).v[1:end-1]...)
    @assert issubset(rng, dom) "$rng SHOULD SUBSETEQ $(dom)"

    # Rename the variables to make an explicit distinction between the input
    # TMs and the loop TMs.
    Dli::Vector{TaylorModelN{numvars,Float64,S}} = Dl0
    Dri::Vector{TaylorModelN{numvars,Float64,S}} = Dr0

    # This rectangle represents the initial set of the current integration
    # iteration. Its time component should always be degenerate: [t, t].
    # FIXME: M. Neher. (2006) p14:
    # "The initial set for the (j + 1)-st integration step is defined by `Rg(Ul,j+1)`".
    # To make explicit that the domain normalization is a form of manual
    # preconditioning, the interval initial set is evaluated after normalization.
    vals::IntervalBox{numvars, S} = IntervalBox([ tm() for tm in Dli ]..., init.v[end])

    println("Dli = "); display(Dli); println()
    println("Dri = "); display(Dri); println()

    # Setup the output vectors / buffers.
    boxes::Vector{IntervalBox} = []
    fboxes::Vector{IntervalBox} = []

    # Round the number of iterations up;
    # the last time step may exceed the time horizon.
    NR_ITERATIONS::Integer = ceil(Int, TIME_HORIZON / TIME_STEP_SIZE)
    # The width threshold below which an interval is considered degenerate.
    # e.g. diam([-1, 1]) = 2 > threshold  =>  NOT degenerate!
    degen_threshold = 1.0e-12

    for it = 1:NR_ITERATIONS
        println("# Start Integration Iteration $it")

        # Get the time component of the initial set.
        tdom = vals.v[end]
        tdiam = diam(tdom)
        @assert tdiam < degen_threshold "The time component of the initial set "*
            "is not considered degenerate: width($tdom) = $tdiam > $degen_threshold."

        # Step 0: Taylorize the dynamics.
        # We want to have a polynomial approximation of the dynamics centered around
        # the midpoint of the current values.
        VarType = typeof(vars[1])
        fpoly = Vector{VarType}(undef, length(vars) - 1)
        vector_field_constructor(fpoly, vcat(vars[1:end-1], (vars[end] + mid(tdom))))
        println("taylorized vector field/dynamics:")
        println(fpoly); println()
        @assert(all([ isassigned(fpoly, idx) for idx in eachindex(fpoly) ]),
            "The dynamics constructor did not assign all vector field components.")

        # Step 1: Obtain the polynomial part of the Taylor model.
	    p::Vector{VarType} = tay_poly(fpoly, k, vars, polynomial.(Dli))
        println("polynomial part of TM:")
        println(p); println()

        # Step 2: Obtain the safe remainder/error interval of the TM.
        #
        # This rectangle is the initial set stretched across the entire time
        # step [ti, ti+δ] of this integration iteration. By definition it only
        # differs from the initial set in its time component.
        #
        # We already shifted t in when Taylorizing the ODEs and computing p.
        # So compensate by setting time to [0, δ]!
        # doms = IntervalBox(vals.v[1:end-1]...,  0..TIME_STEP_SIZE)
        doms = IntervalBox(unitbox(vals.v[1:end-1])..., 0..TIME_STEP_SIZE)

        # FIXME: The real_t parameter can be derived from Dli,
        #        so remove that function argument?
        # Find the safe/contractive remainder.
        safe_rems = tay_model_error(
            vector_field_constructor,
	        p, J, vars, doms, Dli, mid(vals.v[end]),
            NR_CONTRACTIVENESS_TRIES,
            NR_REFINEMENTS,
            REFINEMENT_EPS,
            SCALE)

        # Step 3: Get the new local values (and interval box) and update domain for next step
        # i.e. just change the domain of the time variable in doms
        # Make sure all values are correctly shaped.
        @assert length(p) == length(safe_rems) == length(Dli) == (length(doms)-1)

	# Construct the integrated left Taylor models to get a flowpipe
        Dj = [ TaylorModelN(pj, Ij, doms) for (pj, Ij) in zip(p, safe_rems) ]
        Dri_ext = vcat(Dri, TaylorModelN(vars[end], 0..0, domain(Dri[1])))
	fpipe = IntervalBox([(Dlij(Dri_ext))() for Dlij in Dli]..., tdom.hi..(tdom.hi+TIME_STEP_SIZE))
	vals = IntervalBox([(Dlij(Dri_ext))() for Dlij in Dli]..., interval(tdom.hi+TIME_STEP_SIZE))
	
        # Fix the time variable to the current time; t = ti+δi.
        p = [ pj([vars[1:end-1]..., TaylorN(TIME_STEP_SIZE, k)]) for pj in p ]
        println("Ui as in Neher:")
        println(p)	
        # and construct simpler left Taylor models
        Dj = [ TaylorModelN(pj, Ij, doms) for (pj, Ij) in zip(p, safe_rems) ]
	# to precondition
        Dli, Dri = precondition(Dj, Dri, vars)
        println("Dli = "); display(Dli); println()
        println("Dri = "); display(Dri); println()

        push!(boxes, vals)
        push!(fboxes, fpipe)

        println("                     vals: "); display(vals)
        println("\n=================================\n")
    end

    return boxes, fboxes
end

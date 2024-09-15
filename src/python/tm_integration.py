# Imports
import sympy
import numpy
from sympy.polys.orderings import monomial_key
from math import factorial, floor
from sympy.abc import t, x, y, s

from matplotlib.patches import Rectangle


def create_box(xmin: float, ymin: float, xmax: float, ymax: float) -> Rectangle:
    """Create a rectangle from the given two points.

    The given points (xmin, ymin) and (xmax, ymax) are the lower left and
    upper right corners of the rectangle.
    """
    width = abs(xmax - xmin)
    height = abs(ymax - ymin)
    return Rectangle((xmin, ymin), width, height, linewidth=1, edgecolor='b', fill=False)


def truncate_matrix(matrix: sympy.Matrix, order: int):
  """Truncate all polynomials in the given matrix to the given order.

  This means, for each polynomial all terms of degree > order are dropped.
  Also expands all expressions in the matrix, by distributing multiplications,
  etc.

  Started off adapting an SO post, but ended up mostly consulting the docs.
  https://stackoverflow.com/questions/25001896/sympy-drop-higher-order-terms-in-polynomial
  https://docs.sympy.org/latest/modules/core.html#module-sympy.core.basic
  https://docs.sympy.org/latest/modules/polys/reference.html
  """
  assert order >= 0, "Truncation order must be >= 0."

  matrix = matrix.applyfunc(lambda e: sympy.expand(e))

  truncated = []
  for expr in matrix:
    # All vars in the expression: [x, y, ...]
    vars = list(expr.free_symbols)

    # No variables, meaning the expr is likely 0; continue to prevent exception
    if len(vars) == 0:
      truncated.append(expr)
      continue

    # Dict of the form { term : coefficient }, term is a multivariate monomial
    # e.g. 4x**2 + x*y  =>  { x**2 : 4, x*y : 1 }
    coeffs_dict = expr.as_coefficients_dict()

    # sum(monomial.degree_list) == monomial degree
    poly_trunc = sum(coeff * term for term, coeff in coeffs_dict.items()
                     if sum(sympy.degree_list(term, vars)) <= order)
    truncated.append(poly_trunc)

  return sympy.Matrix(truncated)


def taylor_polynomials(initial_polys: sympy.Matrix,
                       ode_vars: list[sympy.core.symbol.Symbol],
                       ode: sympy.Matrix,
                       order: int) -> sympy.Matrix:
    """Compute the Taylor polynomial over-approximation of an ODE.

    @param[in]  initial_polys  The polynomial to initialize the Lie derivation with.
    @param[in]  oder_vars      The ordered list of ODE variables, e.g. [x, y].
    @param[in]  ode            The ODEs' vector fields as a matrix, e.g.
                               (x'=a, y'=b) becomes [a, b].
    @param[in]  order          The order of the to-generate Taylor polynomials.
    @return     The matrix of generated Taylor polynomials
    """

    # We will be needing the Taylor expansion of the function defined by
    # the dynamics
    lieDers = [initial_polys]  # the 0-th entry is the original TM
    for i in reversed(range(order)):
        J = lieDers[-1].jacobian(sympy.Matrix(ode_vars))
        d = J * ode  # FIXME: missing diff wrt time?
        d = truncate_matrix(d, i)
        lieDers.append(d)

    # We are now ready to compute p1
    taylor_poly = sum([d * (t ** i) / factorial(i) for i, d in enumerate(lieDers)],
            sympy.zeros(*lieDers[0].shape))
    taylor_poly = taylor_poly.applyfunc(
        lambda e: sympy.expand(e, key=monomial_key('lex', ode_vars)))
    return taylor_poly

def picard_operator(taylor_polys: sympy.Matrix,
                    ode: sympy.Matrix,
                    ode_variables: list[sympy.core.symbol.Symbol],
                    order: int) -> sympy.Matrix:
    """Perform the picard operator on the given taylor polynomials pk(x, y, ..., t).

    @param[in]  taylor_polys   The ordered list of taylor polynomials to compute the picard operator of.
    @param[in]  ode            The ODEs w.r.t. which to substitute the polynomials into.
    @param[in]  ode_variables  The ordered list of ODE variables of the given ODEs.
                               e.g. [x, y] for the ODE (x'=..., y'=...).
    @param[in]  order          The Taylor polynomial order to adhere to.
    @return     The matrix of corresponding picard operator results.

    """
    # Substitute polynomial p(x, y, t) into vector field f(x, y).
    substituted = ode.subs(zip(ode_variables, taylor_polys))

    # Integrate the substitution result.
    picard = substituted.applyfunc(lambda e: sympy.integrate(e + s, t))
    picard = picard.applyfunc(lambda e: sympy.expand(e))

    # Perform truncation on the integration result.
    picard = truncate_matrix(picard, order)

    # Prefix the result with the variables to complete the picard operator.
    picard = sympy.Matrix(ode_variables) + picard

    return picard

def refine_remainders(taylor_polys: sympy.Matrix,
                      ode: sympy.Matrix,
                      ode_variables: list[sympy.core.symbol.Symbol],
                      ode_intervals: numpy.ndarray[sympy.plotting.intervalmath.interval],
                      remainder_intervals: numpy.ndarray[sympy.plotting.intervalmath.interval],
                      time: sympy.plotting.intervalmath.interval,
                      order: int) -> numpy.ndarray[sympy.plotting.intervalmath.interval]:
    """Refine a remainder interval vector once, via the picard operator.

    The to-refine remainder intervals do not have to be contractive; this operation
    may also be used to verify contractivenes of the initial remainder estimate vector I0.

    @param[in]  taylor_polys         The ordered list of taylor polynomials to refine the remainders of.
    @param[in]  ode                  The ODEs w.r.t. which to substitute the polynomials into.
    @param[in]  ode_variables        The ordered list of ODE variables of the given ODEs.
                                     e.g. [x, y] for the ODE (x'=..., y'=...).
    @param[in]  ode_intervals        The ordered, initial remainder estimate for each ODE variable.
    @param[in]  remainder_intervals  A running remainder estimate that is refined exactly once during this function.
    @param[in]  time                 The interval valuation of the time variable.
    @param[in]  order                The Taylor polynomial order to adhere to.
    @return     The refined remainders.

    """
    # The time (t) and remainder (s) variables are considered separately from the ODE variables.
    variables = tuple(ode_variables + [t, s])

    # The picard operator facilitates remainder refinement.
    picard = picard_operator(taylor_polys, ode, ode_variables, order)

    # The ode variables take on the same interval value regardless of which component
    # in the taylor polynomial vector.
    interval_values = {
        var: interval
        for var, interval in zip(ode_variables, ode_intervals)
    }

    # The time variable takes on the same interval value regardless of which component
    # in the taylor polynomial vector.
    interval_values.update({ t: time })

    # Apply refinement one-by-one since the s interval valuation can differ by vector component.
    Irefined = []
    for Iprev, exp in zip(remainder_intervals, picard):
        # The (refined) remainder estimate depends on the specific taylor polynomial.
        # i.e. for a vector of taylor polynomials, each may have a different
        # remainder estimate associated with it.
        interval_values.update({ s: Iprev })

        # Perform interval evaluation to compute the actual remainder refinement.
        evaluated_refinement = sympy.lambdify(variables, exp)(
            **{ var.name: interval for var, interval in interval_values.items() }
        )
        Irefined.append(evaluated_refinement)

    return numpy.array(Irefined)

def contractive_remainder(taylor_polys: sympy.Matrix,
                          ode: sympy.Matrix,
                          ode_variables: list[sympy.core.symbol.Symbol],
                          ode_intervals: numpy.ndarray[sympy.plotting.intervalmath.interval],
                          initial_remainders: numpy.ndarray[sympy.plotting.intervalmath.interval],
                          time: sympy.plotting.intervalmath.interval,
                          order: int,
                          nr_contractiveness_tries: int,
                          contractiveness_scale: float,
                          nr_refinements: int) -> numpy.ndarray[sympy.plotting.intervalmath.interval]:
    """Compute a contractive remainder for each pair of true flow and taylor polynomial.

    This operation will first generate a contractive remainder, and then refine the remainder.

    @param[in]  taylor_polys              The ordered list of taylor polynomials to refine the remainders of.
    @param[in]  ode                       The ODEs w.r.t. which to substitute the polynomials into.
    @param[in]  ode_variables             The ordered list of ODE variables of the given ODEs.
                                          e.g. [x, y] for the ODE (x'=..., y'=...).
    @param[in]  ode_intervals             The initial remainder estimate for each ODE variable.
    @param[in]  initial_remainders        The initial remainder estimate I0 that is widened and refined.
    @param[in]  time                      The time slice to find a contractive interval over. Of the form [0, δi]
    @param[in]  order                     The Taylor polynomial order to adhere to.
    @param[in]  nr_contractiveness_tries  The number of times to test an interval for contractiveness.
                                          i.e. the interval will be enlarged at most (#tries - 1) times.
    @param[in]  contractiveness_scale     The factor with which to scale the remainder estimate in case the
                                          contractiveness test fails. e.g. scale = 2.0 and estimate = [-0.1, 0.1],
                                          ==> scale * estimate = [-0.2, 0.2] is tested next
    @param[in]  nr_refinements            The number of times to additionally refine a contractive remainder.
                                          i.e. the contractiveness test maps I0 to I1, meaning one refinement is guaranteed,
                                          after which #refinements additional refinements are applied.
    @return     The contractive, refined remainders.
    """
    assert nr_contractiveness_tries > 0, "Must attempt test contractiveness test at least once."
    assert contractiveness_scale > 1, "Scale must be > 1; subsequent contractiveness tests must enlarge the initial remainder estimate."
    assert nr_refinements >= 0, "Can not apply further refinement a negative number of times."


    uncontractive_remainders = initial_remainders
    Il = None
    for _ in range(nr_contractiveness_tries):
        # I1 is the first "refinement" of the initial/estimate remainder interval vector I0,
        # that correspond to the set of taylor polynomials. The estimate I0 is contractive
        # iff. I1 is a subset of I0 in all of the vector's components.
        I1: numpy.ndarray = refine_remainders(
            taylor_polys, ode, ode_variables, ode_intervals, uncontractive_remainders, time, order)

        contractive: bool = all((
            i0.start <= i1.start <= i1.end <= i0.end
            for i0, i1 in zip(uncontractive_remainders, I1)))

        if contractive:
            Il = I1
            break

        # Contractiveness does not hold, so widen the remainder estimate
        # to hopefully enclose the true flow (the ODE's unique solution).
        uncontractive_remainders = uncontractive_remainders * contractiveness_scale

    if Il is None:
        # The last interval I0 tested for contractiveness did not have the
        # final/redundant widening applied to it, so correct for it in the
        # error message only.
        I1_err = uncontractive_remainders / contractiveness_scale
        raise RuntimeError(f"No contractive remainder found in {nr_contractiveness_tries} "
                            "tries. Last tried contractive test was:"
                           f"\nI0 = {I1_err}"
                           f"\nI1 = {I1}"
                            "\n==> I1[i] subset of I0[i]?")

    for _ in range(nr_refinements):
        # All subsequent refinements may be performed without further contractiveness tests.
        Il = refine_remainders(taylor_polys, ode, ode_variables, ode_intervals, Il, time, order)

    return Il

def initial_set(taylor_polys: sympy.Matrix,
                intervals: numpy.ndarray,
                time_var: sympy.core.symbol.Symbol,
                time: sympy.plotting.intervalmath.interval):
    """Compute the initial set for the next TM integration step.

    @param[in]  taylor_polys  The ordered list of taylor polynomials of the current TM integration interation.
    @param[in]  intervals     The contractive and refined remainder intervals.
    @param[in]  time_var      The sympy variable used as the time variable.
    @param[in]  time          The time step-size δi of the current iteration.
    @return     The initial set for the subsequent TM integration iteration.
    """
    # Given t in [0, δi] = time, fix t=δi in the taylor polynomials.
    # Using evalf minimizes precision errors affecting the result adversely, see docs:
    #   https://docs.sympy.org/latest/modules/core.html#module-sympy.core.evalf
    fixed_time_polys = taylor_polys.evalf(subs={time_var: time.end})

    # The next iteration's initial set is one taylor model per ODE variable:
    # (pk(x, y, ..., t=δi), Il)
    return fixed_time_polys, intervals

def tm_integration_iter(Xl: tuple[sympy.Matrix, list],
                        ode: sympy.Matrix,
                        ode_variables: list[sympy.core.symbol.Symbol],
                        ode_intervals: numpy.ndarray[sympy.plotting.intervalmath.interval],
                        initial_remainders: numpy.ndarray[sympy.plotting.intervalmath.interval],
                        order: int,
                        nr_contractiveness_tries: int,
                        contractiveness_scale: float,
                        nr_refinements: int,
                        time_step_size: float):
    """Perform a single Taylor model integration step.

    @param[in]  Xl                        The input Taylor model initial set for this iteration.
    @param[in]  ode                       The ODEs.
    @param[in]  ode_variables             The ordered list of ODE variables of the given ODEs.
                                          e.g. [x, y] for the ODE (x'=..., y'=...).
    @param[in]  ode_intervals             The initial remainder estimate for each ODE variable.
    @param[in]  initial_remainders        The initial remainder estimate to refine.
    @param[in]  order                     The Taylor polynomial order to adhere to.
    @param[in]  nr_contractiveness_tries  The number of times to test an interval for contractiveness.
    @param[in]  contractiveness_scale     The factor with which to scale the remainder estimate in case the contractiveness test fails.
    @param[in]  nr_refinements            The number of times to additionally refine a contractive remainder.
    @param[in]  time_step_size            The time step-size of the current iteration, i.e. time t in [0, δi] must hold.
    @return     The TM flowpipe generated in this iteration.
    """
    time = sympy.plotting.intervalmath.interval(0.0, time_step_size)

    # STEP 1: Compute the Taylor approximation of the true flow
    taylor_polys = taylor_polynomials(Xl[0], ode_variables, ode, order)

    # STEP 2: Compute a safe remainder estimate.
    Il = contractive_remainder(
        taylor_polys, ode, ode_variables, ode_intervals, initial_remainders,
        time, order, nr_contractiveness_tries, contractiveness_scale, nr_refinements)

    # STEP 3: Compute next iteration's initial set
    Xl = initial_set(taylor_polys, Il, t, time)

    return Xl

def tm_integration(
        ode_variables: list[sympy.core.symbol.Symbol],
        ode: sympy.Matrix,
        ode_intervals: numpy.ndarray[sympy.plotting.intervalmath.interval],
        initial_box: numpy.ndarray[sympy.plotting.intervalmath.interval],
        initial_remainders: numpy.ndarray[sympy.plotting.intervalmath.interval],
        time_horizon: sympy.plotting.intervalmath.interval,
        time_step_size: float,
        time_step_epsilon: float,
        order: int,
        nr_contractiveness_tries: int,
        contractiveness_scale: float,
        nr_refinements: int):
    """Perform TM integration across the entire time_horizon.

    Returns a tuple (flowpipes, boxes), where flowpipes is the ordered
    list of all generated TM flowpipes; and boxes is its associative,
    ordered list of interval enclosures of those TM flowpipes. 

    @param[in]  ode_variables             The ordered list of ODE variables of the given ODEs.
                                          e.g. [x, y] for the ODE (x'=..., y'=...).
    @param[in]  ode                       The ODEs.
    @param[in]  ode_intervals             The initial remainder estimate for each ODE variable.
    @param[in]  initial_box               The domain of all ODE variables at time t=0.
    @param[in]  initial_remainders        The initial remainder estimate to refine.
    @param[in]  time_horizon              The complete time horizon across which to compute flowpipes.
                                          Must be of the form [0, Δ], where Δ is the the time step-sizes
                                          used for each of the flowpipes sum to Δ.
    @param[in]  time_step_size            The time step-size of the current iteration, i.e. time t in [0, δi] must hold.
    @param[in]  time_step_epsilon         An epsilon bound to decide whether the final step-size is large enough to
                                          warrant a flow-pipe. i.e. time horizon [0, 4.1] with step-size 0.2 implies
                                          20 full steps of size [0, 0.2] and a final step of size [0, 0.1]. Choose
                                          epsilon 0.01, then 0.1 > 0.01 meaning we DO want the final step [0, 0.1] to be executed.
                                          If epsilon were 0.2, then 0.1 < 0.2 so that the final step [0, 0.1] would NOT be executed.
                                          This epsilon transparantly handles floating point errors in the time horizon partitioning process.
    @param[in]  order                     The Taylor polynomial order to adhere to.
    @param[in]  nr_contractiveness_tries  The number of times to test an interval for contractiveness.
    @param[in]  contractiveness_scale     The factor with which to scale the remainder estimate in case the contractiveness test fails.
    @param[in]  nr_refinements            The number of times to additionally refine a contractive remainder.
    @return     A tuple (flowpipes, boxes), where flowpipes is the ordered list of all generated TM flowpipes,
                and boxes is its associative, ordered list of interval enclosures of those TM flowpipes. 
    """

    assert len(ode_variables) == len(ode), "Length mismatch: list of free ODE variables does not match number of ODE equations."

    assert time_horizon.start == 0, "Time horizon must start at 0 exactly."
    assert time_horizon.end >= 0, "Time horizon must end >= 0."

    # Determine the number of flowpipes to construct to cover the time horizon.
    nr_full_steps: int = floor(time_horizon.end / time_step_size)
    final_step_size: float = time_horizon.end - nr_full_steps * time_step_size
    # If the final step is too small, then just ignore it.
    # But if the step-size is even smaller, then do it anyways.
    ignore_last_step: bool = final_step_size < time_step_epsilon and\
                             final_step_size < time_step_size
    nr_steps: int = nr_full_steps + (not ignore_last_step)

    #
    # CONSTRUCT FLOWPIPES
    #
    p0 = sympy.Matrix([x, y])   # The polynomial part of the very first initial set X0
    X0 = (p0, initial_remainders)  # The very first initial set X0
    tm_flowpipes = [X0]
    for step in range(nr_steps):
        step_size: float
        if not ignore_last_step and step == nr_full_steps:
            step_size = final_step_size
        else:
            step_size = time_step_size

        # Generate the next flowpipe.
        Xl = tm_flowpipes[-1]
        Xlnew = tm_integration_iter(
            Xl, ode, ode_variables, ode_intervals, initial_remainders,
            order,
            nr_contractiveness_tries,
            contractiveness_scale,
            nr_refinements,
            step_size
        )
        tm_flowpipes.append(Xlnew)


    #
    # CONSTRUCT BOXES
    #

    # The ode variables take on the same interval value regardless of which component
    # in the taylor polynomial vector.
    interval_values = {
        var: interval
        for var, interval in zip(ode_variables, initial_box)
    }

    patches = []
    # Unpacks flowpipe into polynomial pl(Xl, δi) and remainder Il
    for pl, Il in tm_flowpipes:
        # Compute interval enclosure of the flowpipe
        enclosure = sympy.lambdify(ode_variables, pl)(
            **{ var.name: interval for var, interval in interval_values.items() }
        )
        # Peel away redundant layer of nesting: [[ix], [iy], ...] => [ix, iy].
        # Turns a list of one-element lists into just a list of elements.
        if len(enclosure.shape) > 1:
            enclosure = numpy.concatenate(enclosure).ravel()
        
        # The remainder Il represents the uncertainty/error on the polynomial's enclosure.
        enclosure = enclosure + Il

        # Generate 2D box for visualization
        assert enclosure.shape == (2,), "Cannot visualize, expected one 2-dimensional box."
        boxx, boxy = enclosure
        patches.append(create_box(boxx.start, boxy.start, boxx.end, boxy.end))

    return tm_flowpipes, patches

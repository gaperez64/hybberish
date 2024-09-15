# Imports
import sympy
import numpy
from sympy.polys.orderings import monomial_key
from math import factorial, floor
from sympy.abc import t, x, y, s

from matplotlib.patches import Rectangle


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
                       order: int,
                       ode: sympy.Matrix) -> sympy.Matrix:
    """Compute the Taylor polynomial over-approximation of an ODE."""

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
    """Perform the picard operator on the given taylor polynomials pk(x, y, ..., t)."""
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

    initial_remainder: J0, the initial remainder guess that is successively enlarged until
                       contractiveness or failure.
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

    initial_remainder: J0, the initial remainder guess that is successively enlarged until
                       contractiveness or failure.
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
    """Compute the initial set for the next TM integration step."""
    # Given t in [0, δi] = time, fix t=δi in the taylor polynomials.
    # Using evalf minimizes precision errors affecting the result adversely, see docs:
    #   https://docs.sympy.org/latest/modules/core.html#module-sympy.core.evalf
    fixed_time_polys = taylor_polys.evalf(subs={time_var: time.end})

    # The next iteration's initial set is one taylor model per ODE variable:
    # (pk(x, y, ..., t=δi), Il)
    return fixed_time_polys, intervals

def tm_integration_iter(X0: tuple[sympy.Matrix, list],
                        ode: sympy.Matrix,
                        ode_variables: list[sympy.core.symbol.Symbol],
                        ode_intervals: numpy.ndarray[sympy.plotting.intervalmath.interval],
                        initial_remainders: numpy.ndarray[sympy.plotting.intervalmath.interval],
                        order: int,
                        nr_contractiveness_tries: int,
                        contractiveness_scale: float,
                        nr_refinements: int,
                        time_step_size: float):
    """Perform a single Taylor model integration step."""
    time = sympy.plotting.intervalmath.interval(0.0, time_step_size)

    # STEP 1: Compute the Taylor approximation of the true flow
    taylor_polys = taylor_polynomials(X0[0], ode_variables, order, ode)

    # STEP 2: Compute a safe remainder estimate.
    Il = contractive_remainder(
        taylor_polys, ode, ode_variables, ode_intervals, initial_remainders,
        time, order, nr_contractiveness_tries, contractiveness_scale, nr_refinements)

    # STEP 3: Compute next iteration's initial set
    Xl = initial_set(taylor_polys, Il, t, time)

    return Xl

def create_box(xmin: float, ymin: float, xmax: float, ymax: float) -> Rectangle:
    """Create a rectangle from the given two points.
    
    The given points are (xmin, ymin) and (xmax, ymax) are the lower left and
    upper right corners of the rectangle.
    """
    width = abs(xmax - xmin)
    height = abs(ymax - ymin)
    return Rectangle((xmin, ymin), width, height, linewidth=1, edgecolor='b', fill=False)

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

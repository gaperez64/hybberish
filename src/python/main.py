import csv
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from tm_integration import *


if __name__ == "__main__":
    # The initial remainder guesses for ODE variables.
    ode_var_intervals = numpy.array([
        sympy.plotting.intervalmath.interval(-0.0020, 0.0020),  # x  w.r.t. ode_variables
        sympy.plotting.intervalmath.interval(-0.0021, 0.0021)   # y  w.r.t. ode_variables
    ])

    # The initial remainder estimate for the picard operator.
    initial_remainder_estimates = numpy.array([
        sympy.plotting.intervalmath.interval(-0.1, 0.1),  # s variable valuation w.r.t. x in the ode_variables
        sympy.plotting.intervalmath.interval(-0.1, 0.1)   # s variable valuation w.r.t. y in the ode_variables
    ])

    # The initial boxes, i.e. the domain of each ODE variable.
    initial_box = numpy.array([
        sympy.plotting.intervalmath.interval(-1.0, 1.0),  # x  w.r.t. ode_variables
        sympy.plotting.intervalmath.interval(-0.5, 0.5)   # y  w.r.t. ode_variables
    ])

    # The complete time horizon to generate flowpipes over.
    time_horizon = sympy.plotting.intervalmath.interval(0, 0.1)
    # The time step-size of each individual flowpipe
    time_step_size = 0.02
    # When partitioning the time horizon into equal parts, the last
    # part can be of different length compared to all the others.
    # e.g. partition horizon = [0, 1.1] with step-size 0.2 into
    # 5 parts of size [0, 0.2] and one final part of size [0, 0.1].
    # This epsilon is used to decide when that final part is too
    # small to bother with (i.e. the final part is just the result
    # of floating point errors). For example, for a last part [0, 0.001]
    # and epsilon 0.1, then Width([0, 0.001]) = 0.001 < 0.1, so no tm
    # flowpipe is generated for that last section of time horizon.
    # Note: none of the other horizon parts are tested using the epsilon.
    time_step_epsilon = 0.0001

    # Parameters of step 2 of tm integration: remainder refinement
    nr_contractiveness_tries = 10
    contractiveness_scale = 2.0
    nr_refinements = 0

    order: int = 3      # The TM order used for truncation, TM arithmetic, ...
    ode_vars = [x, y]   # Ordered list of all ODE variables
    ode = sympy.Matrix([1 + y, -(x ** 2)])  # The ODEs to over-approximate


    # The ordered list of generated TM flowpipes
    tm_flowpipes: list[tuple[sympy.Matrix, numpy.ndarray]]
    # An associative list of the flowpipes list.
    # The list of interval enclosures of the flowpipes,
    # representing each flowpipe as a box.
    boxes: list

    tm_flowpipes, boxes = tm_integration(
                ode_vars,
                ode,
                ode_var_intervals,

                initial_box,
                initial_remainder_estimates,

                time_horizon,
                time_step_size,
                time_step_epsilon,

                order,
                nr_contractiveness_tries,
                contractiveness_scale,
                nr_refinements)

    coll = PatchCollection(boxes, match_original=True)
    plt.gca().add_collection(coll)
    plt.gca().axis('equal')

    # Output
    plt.savefig('flowpipes.png')

    with open('flowpipes.csv', 'w') as csvfile:
        writer = csv.writer(csvfile, delimiter=',', quotechar='"')
        writer.writerow(("flowpipe idx", "box/enclosure", "polynomial pl", "remainder Il"))

        for i, pair in enumerate(zip(tm_flowpipes, boxes)):
            fp, box = pair
            writer.writerow((i, str(box), None, None))

            polys, intervals = fp
            for poly, interval in zip(polys, intervals):
                writer.writerow((i, None, str(poly), str(interval)))

    plt.show()

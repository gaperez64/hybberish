include("../src/tm_integration_qr.jl")


"""Construct the dynamics.

    This function is used to:
        1. Taylorize the dynamics as input to the Lie derivative routine.
        2. Perform the TM composition in the TM extension of the picard operator.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function f_dot!(du::Vector, u::Vector)
    x, y, t = u
    du[1] = 1.5*x - x*y
    du[2] = -3*y  + x*y
end


# The truncation degree / the degree of all polynomials
# that are used during computations.
ord = 8

# Initial state variable bounds and domain.
initial = [
    ("x", interval(4.9..5.1)),
    ("y", interval(1.9..2.1)),
    ("t", interval(0))
]
init_values = IntervalBox([e[2] for e in initial])

# The fixed time step size.
tstep::Float64 = 0.001
# The number of TM integration algo iterations.
nr_iterations::Integer = 1000
# Specify time as a finite time horizon.
time_horizon::Float64 = nr_iterations * tstep

# The amount of contractiveness attempts before giving up.
nr_contractiveness_tries = 100
# The number of refinements to perform at most.
nr_refinements = 100
# Quit refinement early if the improvement a single refinement
# provides falls below this threshold.
refinement_eps = 0.001
# The scale factor with which to widen the initial safe remainder
# estimate when the contractiveness check fails.
scale_factor = 2.0



boxes::Vector{IntervalBox} = []
fboxes::Vector{IntervalBox} = []
boxes, fboxes = tm_integration_QR(
    f_dot!,
    initial,
    ord,
    (-0.1..0.1),
    time_horizon,
    tstep,
    nr_contractiveness_tries,
    nr_refinements,
    SCALE=scale_factor,
    REFINEMENT_EPS=refinement_eps
)



#
# PLOTTING
#

include("../src/euler.jl")
include("../src/plotting.jl")



"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function ode_euler!(du, u, p, t)
	return f_dot!(du, u)
end

# Start Forward Euler in the middle of the variable domains.
euler_init_state = Vector(mid(init_values))

# Evaluate Forward Euler.
euler_step::Float64 = tstep / 10.0
eseries = euler(ode_euler!, time_horizon, euler_step, euler_init_state)
vnames = get_variable_names()

# Actual plotting
pltND1 = plot_boxes_ND(boxes, vnames, title="Initial sets", legend=false)
pltND2 = plot_boxes_ND(fboxes, vnames, title="Flowpipes", legend=false)

# Add the numerical integration curves to the subplots, for each variable.
for idx in eachindex(pltND1.subplots)
    plot!(pltND1.subplots[idx], eseries[idx], label="Numerically Integrated")
    plot!(pltND2.subplots[idx], eseries[idx], label="Numerically Integrated")
end

# Plot (x, y) pairs, instead of separately plotting (t, x) and (t, y).
# The custom plotting function expects the (x, y) pair to be passed as (y, x).
yx_boxes = [ IntervalBox(box[2], box[1]) for box in fboxes ]
xn, yn = vnames
pltND3 = plot_boxes_ND(yx_boxes, [yn, xn], title="xy-Flowpipes", legend=false)
# Add the numerical integration curve as well.
xy_eseries = [ (ex[2], ey[2]) for (ex, ey) in zip(eseries[1], eseries[2]) ]
plot!(pltND3, xy_eseries, label="Numerically Integrated")

# Compose the plots into a grid.
pltND = plot(pltND1, pltND2)
pltND = plot(pltND, pltND3, layout=(2, 1))


println("Show plot ...")
display(pltND)
println("Press ENTER to continue.")
readline()
println("... done.")

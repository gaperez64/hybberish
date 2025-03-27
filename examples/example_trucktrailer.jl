include("../src/tm_integration.jl")


"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function f_dot!(du::Vector, u::Vector)
	x1, y1, theta0, theta1, v0, dtheta0, = u

    L0 = 0.3375  # truck length
    M0 = 0.1     # truck distance to center of mass
    L1 = 0.3     # trailer length
    M1 = 0.06    # trailer distance to center of mass

	# intermediate variables
	beta01 = theta0 - theta1  # angle between truck and trailer
	v1 = v0 * cos(beta01) + M0 * sin(beta01) * dtheta0  # trailer velocity

	# ODEs + zero(u[.]) for type stability
	du[1] = v1 * cos(theta1) + zero(u[1])
	du[2] = v1 * sin(theta1) + zero(u[2])
	du[3] = dtheta0 + zero(u[3])
	du[4] = (v0 * (1/L1)) * sin(beta01) - (M0 / L1) * cos(beta01) * dtheta0 + zero(u[4])
	du[5] = zero(u[5])  # v0 remains constant
	du[6] = zero(u[6])  # dtheta0 remains constant
end


# The truncation degree / the degree of all polynomials
# that are used during computations.
ord = 4

# Initial state variable bounds and domain.
initial = [
    ("x",   0.00..0.00),
    ("y",   0.25..0.25),
    ("θ0",  1.57..1.57),
    ("θ1",  1.57..1.57),
    ("v0",  0.10..0.30),
    ("dθ0", 0.20..0.40),
    ("t",   0.00..0.00),
]
init_values = IntervalBox([e[2] for e in initial])

# The fixed time step size.
tstep::Float64 = 0.01
# The number of TM integration algo iterations.
nr_iterations::Integer = 5
# Specify time as a finite time horizon.
time_horizon::Float64 = nr_iterations * tstep

# The scale factor for when contractiveness fails.
nr_contractiveness_tries = 10
# The number of refinements to perform at most.
nr_refinements = 10
# Quit refinement early if the improvement a single refinement
# provides falls below this threshold.
refinement_eps = 0.001
# The scale factor with which to widen the initial safe remainder
# estimate when the contractiveness check fails.
scale_factor = 2.0



boxes::Vector{IntervalBox} = []
fboxes::Vector{IntervalBox} = []
boxes, fboxes = tm_integration(
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



function ode_euler!(du, u, p, t)
	f_dot!(du, u)
end

# Start Forward Euler in the middle of the variable domains.
euler_init_state = Vector(mid(init_values))

# Evaluate Forward Euler.
euler_step::Float64 = tstep / 10.0
eseries = euler(ode_euler!, time_horizon, euler_step, euler_init_state)

# FIXME: The indexes specify which variables to plot.
# e.g. To plot only x, choose indexes = [1].
#      To plot plot x, y and v0 choose indexes = [1, 2, 5].
TO_PLOT_ODE_VARIABLES = (1:4)
indexes = vcat(TO_PLOT_ODE_VARIABLES, length(eseries))
# Select the results for parameters v0, dθ0.
vnames = get_variable_names()[indexes]
# Forward Euler.
eseries = eseries[indexes]
# Initial set boxes.
boxes = [IntervalBox(box[indexes]...) for box in boxes]
# Flowpipe enclosure boxes
fboxes = [IntervalBox(box[indexes]...) for box in fboxes]

# Plot the variables against time, as initial sets.
pltND1 = plot_boxes_ND(boxes, vnames, sgtitle="vals", legend=false)
# Plot the variables against time, as flowpipe over-approximations.
pltND2 = plot_boxes_ND(fboxes, vnames, sgtitle="fpipe", legend=false)

# Add the Forward Euler curves to the subplots, for each variable.
for idx in eachindex(pltND1.subplots)
    plot!(pltND1.subplots[idx], eseries[idx], label="Stable ODE Forward Euler")
    plot!(pltND2.subplots[idx], eseries[idx], label="Stable ODE Forward Euler")
end
pltND = plot(pltND1, pltND2)

println("Show plot ...")
display(pltND)
println("Press ENTER to continue.")
readline()
println("... done.")

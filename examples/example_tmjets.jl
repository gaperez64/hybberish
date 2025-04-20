include("../src/tm_integration_qr.jl")


"""Construct the dynamics.

    This function is used to:
        1. Taylorize the dynamics as input to the Lie derivative routine.
        2. Perform the TM composition in the TM extension of the picard operator.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function f_dot!(du::Vector, u::Vector)
    y, t = u
    du[1] = y * cos(y)
end

# The truncation degree / the degree of all polynomials
# that are used during computations.
ord = 20

# Initial state variable bounds and domain.
initial = [
    ("y", interval(0..2)),
    ("t", interval(0))
]
init_values = IntervalBox([e[2] for e in initial])

# The fixed time step size.
tstep::Float64 = 0.01
# The number of TM integration algo iterations.
nr_iterations::Integer = 90
# Specify time as a finite time horizon.
time_horizon::Float64 = nr_iterations * tstep

# The amount of contractiveness attempts before giving up.
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

# Plot the initial set and flowpipe boxes.
pltND1 = plot_boxes_ND(boxes, get_variable_names(), title="Initial Sets (vals)", legend=true)
pltND2 = plot_boxes_ND(fboxes, get_variable_names(), title="Flowpipe Overapprox.", legend=true)

dom = initial[1][2] # Get the domain of y
# Generate traces of the ODEs given various initial conditions of y.
states = [
    [ y, 0.0 ]
    for y in dom.lo:0.05:dom.hi
]
for euler_init_state in states
    eseries = euler(ode_euler!, time_horizon, 0.001, euler_init_state)
    plot!(pltND1, eseries[1])
    plot!(pltND2, eseries[1])
end

pltND = plot(pltND1, pltND2, legend=false, title="Order k=$ord")


println("Show plot ...")
display(pltND)
println("Press ENTER to continue.")
readline()
println("... done.")

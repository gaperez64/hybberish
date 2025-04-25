include("../src/tm_integration_qr.jl")


"""Construct the dynamics.

    This function is used to:
        1. Taylorize the dynamics as input to the Lie derivative routine.
        2. Perform the TM composition in the TM extension of the picard operator.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function f_dot!(du::Vector, u::Vector)
    _u, _v, _t = u
    du[1] = _v
    du[2] = _u^2
end


# The truncation degree / the degree of all polynomials
# that are used during computations.
ord = 8

# Initial state variable bounds and domain.
initial = [
    ("u", interval( 0.95,  1.05)),
    ("v", interval(-1.05, -0.95)),
    ("t", interval( 0.00,  0.00))
]
init_values = IntervalBox([e[2] for e in initial])

# The fixed time step size.
tstep::Float64 = 0.01
# The number of TM integration algo iterations.
nr_iterations::Integer = 100
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

# Actual plotting
pltND1 = plot_boxes_ND(boxes, get_variable_names(), title="Initial sets", legend=true)
pltND2 = plot_boxes_ND(fboxes, get_variable_names(), title="Flowpipes", legend=true)

# Add the Forward Euler curves to the subplots, for each variable.
for idx in eachindex(pltND1.subplots)
    plot!(pltND1.subplots[idx], eseries[idx], label="Numerically Integrated", legend=true)
    plot!(pltND2.subplots[idx], eseries[idx], label="Numerically Integrated", legend=true)
end
pltND = plot(pltND1, pltND2, legend=:outertop)

println("Show plot ...")
display(pltND)
println("Press ENTER to continue.")
readline()
println("... done.")

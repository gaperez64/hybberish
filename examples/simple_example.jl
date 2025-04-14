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
    du[1] = -y - sin(t) + cos(t)
end


# The truncation degree / the degree of all polynomials
# that are used during computations.
ord = 10

# Initial state variable bounds and domain.
# y(0) = [1, 1]
# t(0) = [0, 0]
# This is D_0 specifically in the maths.
initial = [
    ("y", interval(0.95..1.05)),
    ("t", interval(0))
]
init_values = IntervalBox([e[2] for e in initial])

# The fixed time step size.
tstep::Float64 = 0.01
# The number of TM integration algo iterations.
nr_iterations::Integer = 500
# Specify time as a finite time horizon.
time_horizon::Float64 = nr_iterations * tstep

# The scale factor for when contractiveness fails.
nr_contractiveness_tries = 10
# The number of refinements to perform at most.
nr_refinements = 5
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
	y, = u
    # This function implements the stable version of these ODEs.
    # The unstable versions is: "y - sin(t) - cos(t)"
	du[1] =  -y - sin(t) + cos(t)
end

#=
    Note that the solution to the ODEs is "y(t) = cos(t)".
    Given that we know the solution y(t) to the ODEs
        y(t) = cos(t)
        y'(t) = d(cos(t))/dt = - sin(t)
    which we verify against the ODEs
        y'(t) = -y(t) - sin(t) + cos(t)
            = -(cos(t)) - sin(t) + cos(t)
            = - sin(t)
=#

# Start Forward Euler in the middle of the variable domains.
euler_init_state = Vector(mid(init_values))

# Evaluate Forward Euler.
euler_step::Float64 = tstep / 10.0
eseries = euler(ode_euler!, time_horizon, euler_step, euler_init_state)

# Actual plotting
pltND1 = plot_boxes_ND(boxes, get_variable_names(), sgtitle="Initial Sets (vals)", legend=true)
pltND2 = plot_boxes_ND(fboxes, get_variable_names(), sgtitle="Flowpipe Overapprox.", legend=true)
plot!(pltND1, eseries[1], label="Stable ODE Forward Euler")
plot!(pltND1, cos, label="cos(t)") # The ODE solution is "y(t) = cos(t)"
plot!(pltND2, eseries[1], label="Stable ODE Forward Euler")
plot!(pltND2, cos, label="cos(t)") # The ODE solution is "y(t) = cos(t)"
pltND = plot(pltND1, pltND2)

println("Show plot ...")
display(pltND)
println("Press ENTER to continue.")
readline()
println("... done.")


savefig(pltND,  "QR_1_tstep=$(tstep)_#iterations=$nr_iterations.svg")
savefig(pltND1, "QR_2_tstep=$(tstep)_#iterations=$nr_iterations.svg")
savefig(pltND2, "QR_3_tstep=$(tstep)_#iterations=$nr_iterations.svg")

include("../src/tm_integration.jl")



"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function vector_field!(du::Vector, u::Vector)
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

function ttintegration(;
        x::Interval   = 0.00..0.00,
        y::Interval   = 0.25..0.25,
        θ0::Interval  = 1.57..1.57,
        θ1::Interval  = 1.57..1.57,
        v0::Interval  = 0.10..0.30,
        dθ0::Interval = 0.20..0.40,
        t::Interval   = 0.00..0.00,
        truncation_degree::Integer=4,
        time_step_size::Float64=0.01,
        nr_integration_iterations::Integer=1)
    # The truncation degree / the degree of all polynomials
    # that are used during computations.
    ord = truncation_degree

    # Initial state variable bounds and domain.
    initial = [
        ("x",   x),
        ("y",   y),
        ("θ0",  θ0),
        ("θ1",  θ1),
        ("v0",  v0),
        ("dθ0", dθ0),
        ("t",   t),
    ]
    init_values = IntervalBox([e[2] for e in initial])

    # The fixed time step size.
    tstep::Float64 = time_step_size
    # The number of TM integration algo iterations.
    nr_iterations::Integer = nr_integration_iterations
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
    try
        # Suppress the 'println' statements 
        redirect_stdout(devnull)

        boxes, fboxes = tm_integration(
            vector_field!,
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

        redirect_stdout(stdout)

    catch e
        redirect_stdout(stdout)
        throw("TM Integration failed with the following exception: $e")
    end

    return boxes, fboxes
end




# If the script is called from the CLI, then run this code.
if abspath(PROGRAM_FILE) == @__FILE__

    println("Call TM integration ...")
    # TODO: Call this function to compute one trucktrailer initial set!
    boxes::Vector{IntervalBox},
    fboxes::Vector{IntervalBox} = ttintegration()
    println("... done!")

    @assert length(boxes) > 0 "Expected at least one box as output."
    box = boxes[begin]

    println("### The computed initial set ###"); display(box); println()
end


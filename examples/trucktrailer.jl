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

"""Run a number of trucktrailer TM integration iterations.

    @return The vector of generated initial set rectangles.
            Return ::Nothing in case of an exception.
"""
function ttintegration(;
        x::Interval       = 0.00..0.00,
        y::Interval       = 0.25..0.25,
        theta0::Interval  = 1.57..1.57,
        theta1::Interval  = 1.57..1.57,
        v0::Interval      = 0.10..0.30,
        dtheta0::Interval = 0.20..0.40,
        t::Interval       = 0.00..0.00,
        truncation_degree::Integer=4,
        time_step_size::Float64=0.01,
        nr_integration_iterations::Integer=1)
    # The truncation degree / the degree of all polynomials
    # that are used during computations.
    ord = truncation_degree

    # Initial state variable bounds and domain.
    initial = [
        ("x",       x),
        ("y",       y),
        ("theta0",  theta0),
        ("theta1",  theta1),
        ("v0",      v0),
        ("dtheta0", dtheta0),
        ("t",       t),
    ]
    init_values = IntervalBox([e[2] for e in initial])

    # The fixed time step size.
    tstep::Float64 = time_step_size
    # The number of TM integration algo iterations.
    nr_iterations::Integer = nr_integration_iterations
    # Specify time as a finite time horizon.
    time_horizon::Float64 = nr_iterations * tstep

    # The number of times to reattempt the contractiveness test if it fails.
    nr_contractiveness_tries = 10
    # The number of safe remainder refinements to perform at most.
    nr_refinements = 10
    # Quit refinement early if the improvement a single refinement
    # provides falls below this threshold.
    refinement_eps = 0.001
    # The scale factor with which to widen the initial safe remainder
    # estimate when the contractiveness check fails.
    scale_factor = 2.0


    # The list of initial sets produced by TM integration.
    # Each initial set is represented as a box / rectangle.
    initial_sets::Vector{IntervalBox} = []
    try
        # Suppress the 'println' statements 
        redirect_stdout(devnull) do
            initial_sets, _ = tm_integration(
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
        end

    catch e
        @error "TM Integration failed with the following exception: $e"
        return nothing
    end

    return initial_sets
end

"""Run one singular trucktrailer TM integration iteration.

    @return The components (x, y, θ0, θ1) of the result initial set rectangle.
            Return ::Nothing in case of an exception.
"""
function ttintegration(
        x::Interval,
        y::Interval,
        theta0::Interval,
        theta1::Interval,
        v0::Interval,
        dtheta0::Interval,
        t::Interval,
        time_step_size::Float64)
    # Run the more general version of this function to generate the
    # singular initial set rectangle for ONE iteration.
    boxes = ttintegration(
        x=x,y=y,theta0=theta0,theta1=theta1,v0=v0,dtheta0=dtheta0,t=t,
        time_step_size=time_step_size,
        nr_integration_iterations=1)

    if boxes isa Nothing
        return nothing
    end

    # Extract the components (x, y, θ0, θ1) of the first / resulting
    # initial set rectangle.
    return IntervalBox(boxes[1][1:4]...)
end

"""Parse the ARGS builtin global, expecting all arguments to be intervals
    of the form a..b where a and b are floats.

    @param[out] input The parsed input intervals. Should be initialized with
                      default values.
"""
function parse_input_intervals!(input::Vector{I})::Nothing where {I<:Interval}
    args = ARGS

    # There cannot be more input args than the number of variables
    # that tructrailer supports.
    if length(ARGS) > length(input)
        @warn "Truncating redundant input arguments."
        args = args[1:length(input)]
    end

    # Parse string representations of intervals into actual intervals.
    #
    # Strip all whitespace characters.
    args = [replace(a, ' ' => "") for a in ARGS]
    # Parse the inputs as intervals.
    for (idx, a) in enumerate(args)
        bounds = split(a, "..")
        @assert(length(bounds) == 2,
            "Could not split an input on '..' to obtain the lo and hi.")
        lo, hi = bounds
        input[idx] = interval(parse(Float64, lo), parse(Float64, hi))
    end
end




# If the script is called from the CLI, then run this code.
if abspath(PROGRAM_FILE) == @__FILE__

    # Specify default values of the variable domains for the CLI.
    input::Vector{Interval{Float64}} = [
        0.00..0.00, # x
        0.25..0.25, # y
        1.57..1.57, # θ0
        1.57..1.57, # θ1
        0.10..0.30, # v0
        0.20..0.40, # dθ0
        0.00..0.00  # t
    ]

    if length(ARGS) > 7
        @error "Redundant arguments are not allowed, to prevent typos in the "*
               "input. Specify at most 7 interval arguments. Only leave "*
               "whitespace between input arguments."
        exit(0)
    end

    # Parse the CLI arguments and overwrite zero or more elements of the
    # default input values.
    try
        parse_input_intervals!(input)

        # Run TM integration: obtain the initial set for the first four variables.
        result = ttintegration(input..., 0.01)
        @assert result !== nothing "TM integration failed."
        x, y, theta0, theta1 = result
        println("(x, y, θ0, θ1) = $((x, y, theta0, theta1))")
    catch e
        @error "Input parsing OR TM integration failed: $e"
    end
end


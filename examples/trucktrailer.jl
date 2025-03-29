include("../src/tm_integration.jl")
using ArgParse
using CSV

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

"""Run one singular trucktrailer TM integration iteration.

    @return The components (x, y, θ0, θ1, v0, dθ0, t) of the result initial
            set rectangle. Return ::Nothing in case of an exception.
"""
function ttintegration(
        x::Interval,
        y::Interval,
        theta0::Interval,
        theta1::Interval,
        v0::Interval,
        dtheta0::Interval,
        t::Interval;
        truncation_degree::Integer=4,
        time_step_size::Float64=0.01)::IntervalBox
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

    # The fixed time step size.
    tstep::Float64 = time_step_size
    # Specify time as a finite time horizon.
    time_horizon::Float64 = 1 * tstep

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

    # Only output the final initial set, at the end of the time horizon.
    return initial_sets[end]
end

"""Call ttintegration for given input and output files.

    @param[in] ifile The input CSV file path containing points to integrate from.
    @param[in] ofile The output file path to dump the generated points to.
    @param[in] tstep The time step size for TM integration.
    @param[in] rows Process only the first `rows` rows of the input.
                     If rows is ::Nothing, then process the entire input.
"""
function ttintegration(ipath::String, opath::String, tstep::Float64; rows=nothing)::Nothing
    # The input file must exist.
    !isfile(ipath) && return nothing
    # The output dir must exist, but the output file may not yet.
    !isdir(dirname(opath)) && return nothing

    # Create the output file if it does not exist.
    touch(opath)

    # Transform each point.
    reader = CSV.File(ipath)
    ofile = open(opath, "w")

    # Retain the input header; dump the header of the input file to disk.
    header = join(string.(propertynames(reader)), ",") * "\n"
    write(ofile, header)

    # Only process the first `rows` rows of the input.
    # - `rows = 0` is allowed
    # - `rows > length(reader)` defaults to `rows = length(reader)`
    rows = rows === nothing ? length(reader) : min(rows, length(reader))

    for data in reader[1:rows]
        @assert(length(data) == 8,
            "The CSV has an unexpected number of columns: got $(length(data)).")

        # Dump the original data point to disk.
        CSV.write(ofile, [data], writeheader=false, append=true)

        # Unpack the input data to make its column format/ordering explicit.
        x, y, theta0, theta1, v0, dtheta0, v1, t = data

        combinations = []
        try
            # NOTE: Only overwrite variables that change;
            # v0 and dtheta0 remain unchanged after integration.
            x, y, theta0, theta1, _, _, t = ttintegration(
                interval(x), interval(y), interval(theta0), interval(theta1),
                interval(v0), interval(dtheta0), interval(t),
                time_step_size=tstep)

            # Assume time is degenerate.
            t = t.lo
        catch e
            # If TM integration fails, just don't generate additional points.
            continue
        end

        # Generate one border point using each (inf, mid, sup) combination.
        product = Iterators.product # An alias.
        # FIXME: Here we just propagate the input v1. But instead we should
        #        output the v1 computed during TM integration in the vector
        #        field???
        combinations = [
            [ op1(x), op2(y), op3(theta0), op4(theta1), v0, v1, dtheta0, t ]
            for (op1, op2, op3, op4) in product(fill([inf, mid, sup], 4)...)
        ]
        # Transform a vector of rows into a vector of columns.
        combinations = [collect(column) for column in zip(combinations...)]
        # Transform the data to a format that the CSV RowWriter accepts.
        combinations = NamedTuple([
            key => val
            for (key, val) in zip(propertynames(reader), combinations)
        ])

        # Dump the additionally generated data points to disk.
        CSV.write(ofile, combinations, writeheader=false, append=true)
    end
end

"""
    parse_commandline()

    Parse the command-line arguments passed to the program.

    This will return a `Dict` with the following keys:
    - `debug` a `Bool` flag, if true then output debug logging,
"""
function parse_commandline()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "--debug"
            help = "Enable/disable debug logging."
            action = :store_true
        "tstep"
            help = "The fixed time step size for TM integration."
            arg_type = Float64
            required = true
        "input"
            help = "The path to the CSV containing a list of inputs."
            arg_type = String
            required = true
        "output"
            help = "The path to the output CSV to dump the results in."
            arg_type = String
            required = true
        "rows"
            help = "Only process the first `rows` rows of the input."
            arg_type = Int
            default  = nothing
            required = false
    end

    return parse_args(settings)
end



# If the script is called from the CLI, then run this code.
if abspath(PROGRAM_FILE) == @__FILE__

    parsed_args = parse_commandline()
    display(parsed_args)

    # Handle the file input + file output use case.
    ipath_ = parsed_args["input"]
    opath_ = parsed_args["output"]
    tstep_ = parsed_args["tstep"]
    rows_ = parsed_args["rows"]
    ttintegration(ipath_, opath_, tstep_, rows=rows_)
end


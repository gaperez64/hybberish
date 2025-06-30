using ReachabilityAnalysis
using ArgParse
using CSV

@taylorize function vector_field!(du, u, p, t)
    x1, y1, theta0, theta1, v0, dtheta0, = u

    local L0 = 0.3375  # truck length
    local M0 = 0.1     # truck distance to center of mass
    local L1 = 0.3     # trailer length
    local M1 = 0.06    # trailer distance to center of mass

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

function ttintegration(initial, time_horizon)
    # Initial state variable bounds and domain.
    prob = @ivp(x' = vector_field!(x), x(0) ∈ initial, dim=6)

    sol = solve(prob, tspan=(0.0, time_horizon), alg=TMJets21a());


    # Only output the final initial set, at the end of the time horizon.
    # NOTE: The last f() evaluates the vector of taylor models on their domain
    return [f(domain(f)) for f in evaluate(sol[end], tend(sol[end]))]
end

function ttintegration(ipath::String, opath::String)::Nothing
    @assert isfile(ipath) "The input file path does not exist: $ipath"
    @assert !isempty(opath) "The output file path is empty."
    @assert((isdir(dirname(opath)) || isempty(dirname(opath))),
        "The directory '$(dirname(opath))' in output path '$opath' does not exist.")

    thor = 0.5
    @assert thor > 0.0 "Time horizon must have strictly positive value."

    # Create the output file if it does not exist.
    touch(opath)

    # Transform each point.
    reader = CSV.File(ipath)
    ofile = open(opath, "w")

    # Retain the input header; dump the header of the input file to disk.
    header = join(string.(propertynames(reader)), ",") * "\n"
    write(ofile, header)

    rows = length(reader)

    println("\nIntegrating with time horizon $thor")

    for (i, data) in enumerate(reader[1:rows])
        @assert(length(data) == 8,
            "The CSV has an unexpected number of columns: got $(length(data)).")

        println("Processing $i/$rows")

        # Dump the original data point to disk.
        CSV.write(ofile, [data], writeheader=false, append=true)

        # Unpack the input data to make its column format/ordering explicit.
        idx, x, y, theta0, theta1, v0, dtheta0, t = data

        combinations = []
        try
            # NOTE: Only overwrite variables that change;
            # v0 and dtheta0 remain unchanged after integration.
            x0 = Hyperrectangle(low=[x, y, theta0, theta1, v0 - 0.001, dtheta0 - 0.05],
                                high=[x, y, theta0, theta1, v0 + 0.001, dtheta0 + 0.05])
            # Why v0 - 0.01..v0 + 0.01?  # an extra cm/s is fine
            # Why dtheta0 - 0.05..dtheta0 + 0.05?  # an extra 3 deg/s is fine
            res = ttintegration(x0, thor)
            x, y, theta0, theta1, _, _, = res

            # Assume time is degenerate.
            t = t + thor
        catch e
            # If TM integration fails, just don't generate additional points.
            continue
        end

        # Generate one border point using each (inf, mid, sup) combination.
        product = Iterators.product # An alias.
        combinations = [
        [ (idx + 1), op1(x), op2(y), op3(theta0), op4(theta1), v0, dtheta0, t ]
            #f or (op1, op2, op3, op4) in product(fill([inf, mid, sup], 4)...)
        for (op1, op2, op3, op4) in product([inf, sup], [inf, sup], [mid], [mid])
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

function parse_commandline()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "input"
            help = "The path to the CSV containing a list of inputs."
            arg_type = String
            required = true
        "output"
            help = "The path to the output CSV to dump the results in."
            arg_type = String
            required = true
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
    ttintegration(ipath_, opath_)
end


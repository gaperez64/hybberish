using ArgParse, Plots, CSV


"""Generate a plot comparing the input CSV against the generated output CSV.

    @param[in] ifile The input CSV file path containing the to-integrate points.
    @param[in] ofile The output CSV file path containing the integrated points.
    @param[in] irows Only plot the first `irows` rows of the input file.
    @param[in] orows Only plot the first `orows` rows of the output file.
"""
function ttplotting(
        ipath::String,
        opath::String,
        irows::Union{Nothing,Int},
        orows::Union{Nothing,Int})::Nothing
    @assert isfile(ipath) "The input file path does not exist: $ipath"
    @assert isfile(opath) "The output file path does not exist: $opath"

    println("Generating plot ...")

    # Read the input CSV file.
    ifile = CSV.File(ipath)
    ofile = CSV.File(opath)

    # If row truncation is not specified, default to plotting all rows.
    irows = irows isa Nothing ? length(ifile) : min(irows, length(ifile))
    orows = orows isa Nothing ? length(ofile) : min(orows, length(ofile))

    @assert irows >= 0 "Cannot plot a negative number of input rows."
    @assert orows >= 0 "Cannot plot a negative number of output rows."

    # Get the CSV header's column names.
    icolnames = string.(propertynames(ifile))
    ocolnames = string.(propertynames(ofile))
    colnames = icolnames
    @assert(length(icolnames) == length(ocolnames) == INPUT_COLUMN_NR,
        "Detected an invalid number of CSV columns: expected $INPUT_COLUMN_NR")

    # For specific input columns, plot all unique 2D (t, x) points
    # where x is the variable represented by that column and t is time.
    plots = []
    for var_idx in 1:4
        time_idx = INPUT_COLUMN_NR # Assume we the time column's index.
        # Get all unique (t, x) data points from the input/output files.
        # may be duplicate (t, x) points because to create the output file
        # we generate a permutation of > 2 variables.
        unique_points_if = unique([ (r[time_idx], r[var_idx]) for r in ifile[1:irows] ])
        unique_points_of = unique([ (r[time_idx], r[var_idx]) for r in ofile[1:orows] ])

        plt = plot(unique_points_of,
            xlabel=colnames[time_idx],
            ylabel=colnames[var_idx],
            label="CSV points OUT", seriestype=:scatter, markersize=1)
        plot!(plt, unique_points_if, label="CSV points IN", seriestype=:line,
            legendfontsize=8, background_color=:transparent, legend=:outertop)

        push!(plots, plt)
    end
    # Compose the subplots.
    pltND = plot(plots..., background_color=:white)

    # Save the plot by adapting the output CSV path.
    savefig(pltND, opath * ".svg")

    println("... done!")
end

"""
    parse_commandline()

    Parse the command-line arguments passed to the program.

    This will return a `Dict`.
"""
function parse_commandline()
    settings = ArgParseSettings()

    @add_arg_table! settings begin
        "input"
            help = "The path to the input CSV containing input data points."
            arg_type = String
            required = true
        "output"
            help = "The path to the output CSV containing integrated points."
            arg_type = String
            required = true
        "--irows"
            help = "Only process the first `irows` rows of the input."
            arg_type = Int
            default  = nothing
            required = false
        "--orows"
            help = "Only process the first `orows` rows of the output."
            arg_type = Int
            default  = nothing
            required = false
    end

    return parse_args(settings)
end


# We expact an exact number of input columns.
# Enforcing this requirement ensures the input CSV is generally formatted
# as we expect.
INPUT_COLUMN_NR::Int = 7

# If the script is called from the CLI, then run this code.
if abspath(PROGRAM_FILE) == @__FILE__

    parsed_args = parse_commandline()
    display(parsed_args)

    # Handle the file input + file output use case.
    ipath_ = parsed_args["input"]
    opath_ = parsed_args["output"]
    irows_ = parsed_args["irows"]
    orows_ = parsed_args["orows"]

    ttplotting(ipath_, opath_, irows_, orows_)
end

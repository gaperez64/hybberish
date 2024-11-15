include("plotting.jl")
include("euler.jl")
include("scripts/clean_intervals.jl")
include("tm_integration.jl")

if abspath(PROGRAM_FILE) == @__FILE__

    # Full TM integration, based on Example 3.3.2 and Example 3.3.6
    k = 3   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 5
    NR_REFINEMENTS           = 1
    SCALE                    = 2.0
    TIME_STEP_SIZE = 0.02
    TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
    USE_LOCAL_HORIZON = true
    time_horizon = 0.8     # The finite time horizon
    vars = set_variables("x y t", order=k)
    vars_no_t = get_variable_names()[1:end-1]
    # The vector field f of the ODEs:
    #   f[1] = 1 + y
    #   f[2] = -x^2
    f = [1 + vars[2],  # x
        -vars[1]^2]   # y
    dom = IntervalBox([-1..1,      # x
                -0.5..0.5,  # y
                0..0.02])   # t
    # Initial remainder estimate J, a hyperrectangle
    J = fill(-0.1..0.1, length(f))


    initial_sets, step_sizes =
        tm_integration(f, dom, k, J, time_horizon,
                       TIME_STEP_SIZE,
                       TIME_STEP_SIZE_EPS,
                       NR_CONTRACTIVENESS_TRIES,
                       NR_REFINEMENTS,
                       SCALE)
    println("\ninit sets ($(length(initial_sets))):")
    for e in initial_sets
        println("    $e")
    end
    println("\nδi vector ($(length(step_sizes))) = $step_sizes")

	initial_boxes = map((Xi) -> interval_initial_set(Xi, dom), initial_sets)
	# FIXME: Assume t is the last variable, and drop its dummy vector element.
	initial_boxes = map((Bij) -> Bij[1:end-1], initial_boxes)

    println("\ninit boxes ($(length(initial_boxes))):")
    for e in initial_boxes
        println("    $e")
    end

    plt3 = plot_boxes_ND(initial_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                         title="Computed boxes Bi", titlefontsize=8)

    # Create a sequence of known boxes that have some overlap with
    # the initial boxes.
    known_boxes = []
    for box in initial_boxes
        known_box = []
        for interval in box
            known_int = interval - diam(interval) / 2.0
            push!(known_box, known_int)
        end
        push!(known_boxes, known_box)
    end
    known_boxes = get_clean_boxes() # FIXME: insert the clean boxes
    plt4 = plot_boxes_ND(known_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                         title="Given boxes Gi", titlefontsize=8)

    plt5 = nothing
    if USE_LOCAL_HORIZON
        # Use a hack to plot TWO sequences of boxes on the same figure:
        # Append the box vectors to each other. BUT, a shifted version of the
        # time step vector, since the global time in the plotting function
        # increases for each box plotted.
        plt5 = plot_boxes_ND([initial_boxes..., known_boxes...], [step_sizes..., step_sizes...], vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                            title="Bi & Gi")
    else
        plt5 = plot_boxes_ND([initial_boxes..., reverse(known_boxes)...], [step_sizes..., (-step_sizes)...], vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                            title="Bi & Gi")
    end

    intersect_boxes = map(((init, known),) -> intersect(IntervalBox(init), IntervalBox(known)), zip(initial_boxes, known_boxes))
    intersect_not_empty = map((int) -> emptyinterval() != int, intersect_boxes)


    plt6 = plot_boxes_ND(intersect_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                         title="Bi ∩ Gi")

    figures = [plt3, plt4, plt5, plt6]
    # Find the single ylim to enclose all plots.
    ylims_val = hull(Interval.(ylims.(figures))...)
    ylims_val = (ylims_val.lo, ylims_val.hi)
    xlims_val = hull(Interval.(xlims.(figures))...)

    # We want there to be overlap between the computed and known boxes.
    @assert all(intersect_not_empty)

    plt_composed = plot(figures..., layout=(1, length(figures)),
                        ylims=ylims_val, xticks=[xlims_val.lo, mid(xlims_val), xlims_val.hi])

    euler_step_size = 0.001
    datax_max = time_horizon
    tseries, yseries = euler(ode!, datax_max, euler_step_size, [0.0; 0.0])
    println("euler t: $tseries")
    println("euler y: $yseries")
    datax = USE_LOCAL_HORIZON ? [0.0, map((_) -> euler_step_size, tseries[2:end])...] : tseries
    datay = map((vec) -> vec[2], yseries)

    plot_recursively!(plt_composed, datax, [datay, datay], vars_no_t)

    png_dst = "output/Xi_composed_local.png"
    println("Press Enter to continue...")
    display(plt_composed)
    readline()
    savefig(plt_composed, png_dst)
    println("Wrote .png to \"$png_dst\"")

    plt_2D = plot_boxes_2D(initial_boxes, vars_no_t)
    datax = map((vec) -> vec[1], yseries)
    datay = map((vec) -> vec[2], yseries)

    plot_recursively!(plt_2D, datax, datay)

    
    png_dst = "output/Xi_2d_local.png"
    println("Press Enter to continue...")
    display(plt_2D)
    readline()
    savefig(plt_2D, png_dst)
    println("Wrote .png to \"$png_dst\"")
end

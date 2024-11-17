using Plots
using TaylorModels # Includes IntervalArithmetic


"""Construct a rectangle shape for plotting.

    @param[in] w The width of the rectangle
    @param[in] h The height of the rectangle
    @param[in] x The x-coord of the bottom left corner of the rectangle
    @param[in] y The y-coord of the bottom left corner of the rectangle
"""
rect(w, h, x, y) = Shape(x .+ [0, w, w, 0, 0], y .+ [0, 0, h, h, 0])


"""For a sequence of 2-dimensional boxes, plot the box's components against each other.

    For a sequence of boxes [(Ix, Iy), ...], plot interval Ix against
    interval Iy for each 2D box (Ix, Iy).

    @param[in] boxes          The sequence of boxes to plot.
    @param[in] variable_names The ordered name of the ODE variables.
    @param[in] title             The plot title.
    @return The generated plot object.
"""
function plot_boxes_2D(boxes, variable_names::Array{String}; title::String="", titlefontsize::Integer=10)
    # Require there to be exactly two variables.
    @assert length(variable_names) == 2

    xname, yname = variable_names
    rectangles = map((box) -> rect(diam(box[1]), diam(box[2]), box[1].lo, box[2].lo), boxes)
    return plot(rectangles, fc=:transparent, lc=:blue, legend=:false, xlabel=xname, ylabel=yname, title=title, titlefontsize=titlefontsize)
end


"""For a sequence of n-dimensional boxes, plot all components of each box against
    time individually.

    For a sequence of boxes [(I1, ..., In), ...], plot intervals I1, ..., In
    against time individually, for each nD box (I1, ..., In).

    @param[in] boxes             The sequence of boxes to plot.
    @param[in] time_steps        The time steps δ1, ..., δN. The time
                                 step δi corresponds to the i-th box.
    @param[in] variable_names    The ordered name of the ODE variables.
    @param[in] use_local_horizon If true, then plot each box against
                                 `t in [0, δi]`. Else, for
                                 `T = summ_{j=1}^{i-1} δj` plot each box
                                 against `t in [T, T + δi]`
    @param[in] title             The plot title.
    @return The generated plot object.
"""
function plot_boxes_ND(boxes, time_steps::Array{Float64}, variable_names::Array{String};
                       use_local_horizon::Bool=false, title::String="", titlefontsize::Integer=10)
    # Each box must correspond to a time step.
    @assert length(boxes) == length(time_steps)

    nr_ODE_components::Integer = length(variable_names)

    # Compute one plot for each ODE component.
    rectangles = map((_) -> Vector(), variable_names)
    # The time interval to plot each variable against.
    time_step = Interval(0.0, 0.0)
    # For each box, and the corresponding time step δi,
    # generate the corresponding rectangles.
    for (box, step) in zip(boxes, time_steps)
        if use_local_horizon
            time_step = Interval(0.0, step)
        else
            time_step = Interval(time_step.hi, time_step.hi + step)
        end

        # Each ODE component corresponds to a fixed box component
        # and has a corresponding set of shapes (rectangles) to draw.
        for (dimension, shape_set) in zip(box, rectangles)
            shape = rect(time_step.hi - time_step.lo, diam(dimension), time_step.lo, dimension.lo)
            push!(shape_set, shape)
        end
    end

    plots = []
    # For each variable, plot the corresponding rectangles in its own figure.
    for (var, shape_set) in zip(variable_names, rectangles)
        shape_set = map((e) -> e, shape_set)    # Convert Array to Vector, because plotting does not want arrays :(
        plt = plot(shape_set, fc=:transparent, lc=:blue, legend=:false,
                   xlabel="t (time)", ylabel="$var")
        push!(plots, plt)
    end

    # Compose the separate figures into a column of figures.
    return plot(plots..., layout=(nr_ODE_components, 1), title=title, titlefontsize=titlefontsize)
end

"""A flattened list of all leaf subplots of the given plot.

    This is a simple wrapper around the index accessing of
    a plot's subplots, to retrieve all of them as an ordered
    list in the same order as the indexing would access them.
"""
function subplots(plt::Plots.Plot)
    return map((subplot_idx) -> plt[subplot_idx], 1:1:length(plt))
end

"""Plot the given x, y data on the leaf subplots of the specified plot.

    The leaf subplots are modified in-place. Given
               y = [s1, ..., sn]
        y_labels = [l1, ..., ln]
    where si is the i-th y-axis series and li is the corresponding label.
    For any subplot, the first series sj whose label lj matches the y-axis
    label is drawn on that subplot.

    @param[in,out] plt  The plot to modify
    @param[in] x        The x data to plot, a single series
    @param[in] y        The y data to plot, a list of series
    @param[in] y_labels The y-axis labels. The i-th label corresponds to
                        the i-th y series.
"""
function plot_recursively!(plt, x, y, y_labels)
    # Each series in y corresponds to one label.
    @assert length(y) == length(y_labels)

    for subplt in plt.subplots
        for (series, label) in zip(y, y_labels)
            ylabel = subplt[:yaxis][:guide]
            if label == ylabel
                plot!(subplt, x, series, linewidth=2, linecolor = :red)
                break
            elseif label == y_labels[end]
                # This axis' ylabel does not correspond to any series
                @assert(false)
            end
        end
    end
end

"""Plot the given x, y data on the leaf subplots of the specified plot.

    The leaf subplots are modified in-place. The y-axis data
    is drawn to all subplots.

    @param[in,out] plt  The plot to modify
    @param[in] x        The x data to plot, a single series
    @param[in] y        The y data to plot, a single series
"""
function plot_recursively!(plt, x, y)
    for subplt in plt.subplots
        plot!(subplt, x, y, linewidth=2, linecolor = :red)
    end
end

"""Plot the computed boxes, known boxes, their intersection,
    and the approximation of the ODE computed with forward Euler.

    @param[in] initial_boxes The boxes of the computed initial sets
    @param[in] known_boxes The boxes to compare our solution to.
                           e.g. boxes gotten through flowstar
    @param[in] step_sizes The step sizes (deltas) for all boxes
    @param[in] vars_no_t The ordered string names of all variables except time
    @param[in] USE_LOCAL_HORIZON Whether to use the local time horizon for plotting
    @param[in] euler_ode The ODE function to pass to the forward euler routine
    @param[in] euler_init_state The initial conditions to pass to the forward euler routine
"""
function plot_vars_against_time(initial_boxes::Vector{IntervalBox{N, Float64}}, known_boxes::Vector{IntervalBox{N, Float64}},
    time_horizon::Float64, step_sizes::Vector{Float64}, vars_no_t::Vector{String}, USE_LOCAL_HORIZON::Bool,
    euler_ode::Function, euler_init_state::Vector{Float64}) where N
    plt_computed = plot_boxes_ND(initial_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
    title="Computed boxes Bi", titlefontsize=8)

    plt_known = plot_boxes_ND(known_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
    title="Given boxes Gi", titlefontsize=8)

    plt_both = nothing
    if USE_LOCAL_HORIZON
        # Use a hack to plot TWO sequences of boxes on the same figure:
        # Append the box vectors to each other. BUT, a shifted version of the
        # time step vector, since the global time in the plotting function
        # increases for each box plotted.
        plt_both = plot_boxes_ND([initial_boxes..., known_boxes...], [step_sizes..., step_sizes...], vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
        title="Bi & Gi")
    else
        plt_both = plot_boxes_ND([initial_boxes..., reverse(known_boxes)...], [step_sizes..., (-step_sizes)...], vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
        title="Bi & Gi")
    end

    intersect_boxes = map(((init, known),) -> intersect(IntervalBox(init), IntervalBox(known)), zip(initial_boxes, known_boxes))
    intersect_not_empty = map((int) -> emptyinterval() != int, intersect_boxes)

    # We want there to be overlap between the computed and known boxes.
    @assert all(intersect_not_empty)

    plt_intersect = plot_boxes_ND(intersect_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
    title="Bi ∩ Gi")

    figures = [plt_computed, plt_known, plt_both, plt_intersect]

    # Find the single ylim to enclose all plots.
    ylims_val = hull(Interval.(ylims.(figures))...)
    ylims_val = (ylims_val.lo, ylims_val.hi)
    xlims_val = hull(Interval.(xlims.(figures))...)

    plt_composed = plot(figures..., layout=(1, length(figures)),
    ylims=ylims_val, xticks=[xlims_val.lo, mid(xlims_val), xlims_val.hi])

    # Compute forward Euler
    euler_step_size = 0.001
    datax_max = time_horizon
    tseries, vseries = euler(euler_ode, datax_max, euler_step_size, euler_init_state)

    datax = USE_LOCAL_HORIZON ? [0.0, map((_) -> euler_step_size, tseries[2:end])...] : tseries


    @assert length(vseries) > 0 "Forward Euler did not produce any values."
    datay = []
    for subseries_idx in eachindex(vseries[1])
    subseries = map((vec) -> vec[subseries_idx], vseries)
    push!(datay, subseries)
    end

    plot_recursively!(plt_composed, datax, datay, vars_no_t)

    return plt_composed, tseries, vseries
end

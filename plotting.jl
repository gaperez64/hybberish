using Plots
using IntervalArithmetic


"""For a sequence of 2-dimensional boxes, plot the box's components against each other.

    For a sequence of boxes [(Ix, Iy), ...], plot interval Ix against
    interval Iy for each 2D box (Ix, Iy).

    @param[in] boxes          The sequence of boxes to plot.
    @param[in] variable_names The ordered name of the ODE variables.
    @return The generated plot object.
"""
function plot_boxes_2D(boxes::Vector{T}, variable_names::Tuple{String, String};
        kwargs...) where T <: IntervalBox
    # Require there to be exactly two variables.
    @assert all(@. length(boxes) == 2) "All boxes must be 2D."
    @assert !(:xlabel in keys(kwargs)) "Invalid kwarg: $(:xlabel) must" *
        " be passed as a positional param."
    @assert !(:ylabel in keys(kwargs)) "Invalid kwarg: $(:ylabel) must" *
        " be passed as a positional param."

    xname, yname = variable_names
    return plot(boxes; xlabel=xname, ylabel=yname,
        fillalpha=0.0, lc=:blue, legend=:false, kwargs...)
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
    @return The generated plot object.
"""
function plot_boxes_ND(boxes::Vector{T}, time_steps::Array{Float64},
        variable_names::Array{String};
        use_local_horizon::Bool=false, kwargs...) where T <: IntervalBox
    @assert length(boxes) == length(time_steps) "Each box must correspond"*
        " to a time step."

    VAR_TIME_NAME::String = "time"
    time_intervals::Array{Interval} = []

    # Map each time step δi to interval [0, δi].
    if use_local_horizon
        time_intervals = map((di) -> 0..di, time_steps)
    # Map each time step di to interval [T, T + δi].
    else
        cs = cumsum(time_steps)
        cs = zip(vcat(0, cs[1:end-1]), cs)
        time_intervals = Interval.(cs)
    end

    plots = []
    for (vidx, var) in enumerate(variable_names)
        boxes2d::Vector{IntervalBox{2, Float64}} = [
            IntervalBox(time_intervals[bidx], box[vidx])
            for (bidx, box) in enumerate(boxes)
        ]

        push!(plots, plot_boxes_2D(boxes2d, (VAR_TIME_NAME, var)))
    end

    # Compose the separate figures into a column of figures.
    return plot(plots..., layout=(length(variable_names), 1), kwargs...)
end

"""Plot a sequence of n-dimensional (nD) boxes as 2-dimensional (2D) boxes.

    Assume the last component of each box is the time component,
    which is used as the x-axis component. For a sequence of nD boxes
        [(I_1, ..., I_n), ...]
    plot the nD box (I_1, ..., I_n) as n-1 boxes
        (I_n, I_1), ..., (I_n, I_(n-1)).

    @param[in] boxes          The sequence of nD boxes to plot.
    @param[in] variable_names The ordered name of the ODE variables.
    @param[in] kwargs         Kwargs are passed to the generated plot.
    @return The generated plot object.
"""
function plot_boxes_ND(
        boxes::Vector{T},
        variable_names::Array{String};
        kwargs...) where T <: IntervalBox

    @assert all(length.(boxes) .> 1) "Cannot generate 2D plots for 1D boxes."

    nr_vars_no_t = length(variable_names) - 1
    plots = []
    for vidx in 1:nr_vars_no_t
        # Assume the last component of the boxes and variable names
        # correspond to the time dimension / component.
        tname = variable_names[end]
        vname = variable_names[vidx]

        # For the i-th var, plot the i-th box component against time.
        boxes2d = [ IntervalBox(box.v[end], box[vidx]) for box in boxes ]
        plt = plot(
            boxes2d;
            xlabel=tname, ylabel=vname,
            fillalpha=0.0, lc=:blue, legend=:false
        )
        push!(plots,  plt)
    end

    # Compose the separate figures into a column of figures.
    return plot(plots..., layout=(nr_vars_no_t, 1); kwargs...)
end

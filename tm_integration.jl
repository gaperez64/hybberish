using Plots
using TaylorModels
include("flowstar.jl")
include("plotting.jl")
include("euler.jl")
include("scripts/clean_intervals.jl")

"""Convert the given box to a TM initial set.

    A box also represents a TM initial set. This is a different,
    but equivalent representation to a TM initial set.
"""
function tm_initial_set(box::IntervalBox)
    # Each variable should have a matching box element.
    @assert(get_numvars() == length(box))

    # Do the conversion: x  =>  (x, [0, 0])  given the domain of x.
    zd = zero(box)
    vars = get_variables()
    # FIXME: This also adds an entry for variable t, but it should not?
    #       Because an ODE with components x,y should have a 2-dimensional
    #       init set Xi, with components only for x,y and not t?
    return map((var) -> TaylorModelN(var, 0..0, zd, box), vars)
end


"""Convert the given TM initial Xi set to a box Bi.

    A box also represents a TM initial set. This is a different,
    but equivalent representation to a TM initial set.
"""
function interval_initial_set(Xi, domain::IntervalBox)
    # Each variable except t should have a matching TM initial set element.
    @assert(get_numvars() == length(Xi))

    # Do the conversion: (p, I)  =>  p(domains) + I  given the domains.
    Bij(Xij) = polynomial(Xij)(domain) + remainder(Xij)
    return map((Xij) -> Bij(Xij), Xi)
end


"""Given a single TM (p, I), compute the corresponding flowpipe.

    Given an initial set Xi, the flowpipe is computed as (p(Xi, t), I).
"""
function flowpipe(p::TaylorN, I::Interval, Xi::Vector, domain::IntervalBox)
    Xi = copy(Xi)
    # Compute the flowpipe Fi = (p(Xi), I) = (ps, Is + I)
    #             where p(Xi) = (ps, Is)
    tm = p(Xi)
    ps = polynomial(tm)
    Is = remainder(tm)
    zd = zero(domain)
    return TaylorModelN(ps, Is + I, zd, domain)
end


"""Compute the next iteration's initial set Xi based on the
    current flowpipe Fi.

    The time step-size δi is used for substitution in a TaylorSeries.
"""
function initial_set(Fi::TaylorModelN, delta_i::TaylorN, variables, domain::IntervalBox)
    # To compute Xi, fix t=δi in Fi.
    # FIXME: Assume t is the last variable and exclude it from the list.
    vars_except_t = variables[1:end-1]
    valuations = vcat(vars_except_t, [delta_i])
    ps = polynomial(Fi)(valuations)
    Is = remainder(Fi)
    zd = zero(domain)
    return TaylorModelN(ps, Is, zd, domain)
end

"""The standard TM integration algorithm.

    Compute an overapproximation of the true flow of the
    system of ODEs represented by the given vector field,
    as a sequence of flowpipes over partial time horizons
    [0, δi] of the full time horizon [0, Δ].
"""
function tm_integration(f, domain, k::Integer, J, time_horizon::Float64,
                        TIME_STEP_SIZE::Float64,
                        TIME_STEP_SIZE_EPS::Float64,
                        NR_CONTRACTIVENESS_TRIES::Integer,
                        NR_REFINEMENTS::Integer,
                        SCALE::Float64)
    @assert(time_horizon > 0)  # The time horizon must not be [0, 0].
    @assert(TIME_STEP_SIZE > 0)
    @assert(NR_CONTRACTIVENESS_TRIES >= 0)
    @assert(NR_REFINEMENTS >= 0)

    # Derive the corresponding TM initial set from the interval initial set.
    X0 = tm_initial_set(domain)
    zd = zero(domain)

    # Start the TM integration loop.
    remaining_time::Float64 = time_horizon
    initial_sets::Array = [X0]
    step_sizes::Array = [0.0]
    Fi = nothing    # The current flowpipe.
    Xi = X0         # The current initial set.
    # FIXME: Would a for-loop be cleaner?
    while remaining_time > 0.0
        delta_i = min(remaining_time, TIME_STEP_SIZE)
        remaining_time -= delta_i
        println("time step [0, $delta_i] (remaining Δ: $remaining_time)")

        if delta_i < TIME_STEP_SIZE_EPS
            println("==> Skipping! δi = $delta_i < $TIME_STEP_SIZE_EPS, the step-size is too small.")
            continue
        end

        # Step 1: generate the poly approximation of the flow.
        p = tay_poly(f, k)

        # Step2: generate a safe, refined remainder interval.
        I = tay_model_error(f, p, domain, k, J,
                            NR_CONTRACTIVENESS_TRIES,
                            NR_REFINEMENTS,
                            SCALE)

        # Step 3: Compute the flowpipe and next initial set.
        # get a copy of the t variable before messing up order
        t = get_variables()[end]
        old_order = get_order()
        old_variables = get_variables()
        # FIXME: AFAIK the 'evaluate' function for TaylorN of TaylorSeries
        # acts up when the substitution values have different orders.
        # For example:
        #   x, y, a, b = set_variables("x y a b", order=3)
        #   p = TaylorN(x + y)
        #   p([1, 2, 3, 4])   # This works fine.
        #   p([a, b, x, y])   # This works fine too.
        #   p([1, b, 3, y])   # ERROR.
        # This even fails if one of the variables is of lower order than others!
        #   x, y, a, b = set_variables("x y a b", order=3)
        #   p = TaylorN(x + y)
        #   p([a, b, x, y])   # This works fine.
        #   a = get_variables(1)[3]
        #   p([a, b, x, y])   # ERROR. a is now order 1, which is less than 3!
        #
        # FIXME: The solution: Convert the scalar δi to a TaylorN.
        # Construct a TaylorN representation of the scalar δi.
        # Must be of the same order as the variables, else evaluation fails!
        delta_i_tayn = TaylorN(delta_i, get_order())

        # To account for composition, we'll up the current order
        hgr_order = old_order * maximum(get_order.(p))
        vars = set_variables(get_variable_string(), order=hgr_order)

        # TODO: Compute the flowpipe components.
        # FIXME: Does TaylorModels actually add the enclosures of truncated
        # terms during TM arithmetic to the remainders? Because the
        # substitution below requires this!
        # Also, does TaylorModels take into account that TaylorSeries seems
        # to compute the order of a constructed expression based on the order
        # of the lowest order term.
        # e.g. for x order 1 and y order 3, then (x + y^2) is order 1 and
        # actually becomes x with y^2 truncated.
        tmv = zip(p, I)
        Fi = map(((pj, Ij),) -> flowpipe(pj, Ij, Xi, domain), tmv)

        # TODO: Compute the initial set vector.
        Xi::Vector{TaylorModelN{3, Float64, Float64}} =
             map((Fij) -> initial_set(Fij, delta_i_tayn, old_variables, domain), Fi)

        dummy_t_tm = TaylorModelN(t, 0..0, zd, domain)
        push!(Xi, dummy_t_tm)   # FIXME: Append dummy
        # FIXME: output Xi AFTER appending the dummy TM, for consistency
        # with X0 that was added to the output list before the loop.
        push!(initial_sets, Xi)
        push!(step_sizes, delta_i)

        # We now restore the order.
        vars = set_variables(get_variable_string(), order=old_order)
    end

    return initial_sets, step_sizes
end


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

0
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
    known_boxes = CLEAN_BOXES # FIXME: insert the clean boxes
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

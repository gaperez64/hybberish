### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ abf81117-2365-4a81-8116-73baf6f8caaf
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ 2580ed6a-98dc-427b-a5d8-344c5a9ed4ad
begin
	using Plots
	using TaylorModels
	include("plotting.jl")
	include("euler.jl")
	include("scripts/clean_intervals.jl")
	include("tm_integration.jl")
end

# ╔═╡ 6b2d0d14-a438-11ef-11d9-e9ba654fb47e
md"""
# TM integration

This notebook demonstrates the TM integration algorithm with hardcoded inputs.
You can play around with these inputs in the cells below, 

It is based on **Flowstar example 3.3.10**.
"""

# ╔═╡ 45cae1dc-ad20-4a5c-94ec-188b2e0380f2
md"""
## Environment

We first set up the environment, by including the necessary packages and files, etc.
"""

# ╔═╡ 774572a0-cfbe-4d4d-926e-c0964715608e
md"""
## Configure TM integration

We now set up the hardcoded inputs. The comments in the markdown cells will be kept short. For more detailed explanations, see the notebook for example 3.3.6.
"""

# ╔═╡ e9eaca80-5d03-45d0-ace0-9e9939200aa6
begin
    # Full TM integration, based on Example 3.3.2 and Example 3.3.6
    k = 3   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 5
    NR_REFINEMENTS           = 1
    SCALE                    = 2.0
    TIME_STEP_SIZE = 0.02
    TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
    time_horizon = 0.8     # The finite time horizon
    vars = set_variables("x t", order=k)
    vars_no_t = get_variable_names()[1:end-1]

    dom = IntervalBox([-1..1,      # x
    0..0.02])   # t
    vars_tms = tm_initial_set(dom)
    
    
    # FIXME: Bumping up the order to avoid assertion errors
    # during TM arithmetic.
    old_order = get_order()
    new_order = old_order*2
    vars = set_variables(get_variable_string(), order=new_order)
    
    # The vector field f of the ODEs:
    #   f[1] = sin(x)
    # Represent the vector field as TMs:
    #   f[1] = x - (1/6)x^3 + I
    # where I is the Lagrange remainder of approximating
    # `sin(x)` by `x - (1/6)x^3`
    f = [ sin(vars_tms[1]) ]   # x
    # Initial remainder estimate J, a hyperrectangle
    J = fill(-0.1..0.1, length(f))

    # TODO: Reset to old order.
    vars = set_variables(get_variable_string(), order=old_order)

	f # Pluto cell output
end

# ╔═╡ 73c9ed6c-391c-4628-809b-60ec755afb94
md"""
## Run TM integration

iteration |   -   |      1     |   $\dots$  | N
    :---  | :---: |    :---:   |     :---:  | :---:
init set  | $X_0$ |    $X_1$   |   $\dots$  | $X_N$
step size |   0   | $\delta_1$ |   $\dots$  | $\delta_N$
"""

# ╔═╡ 7ec7bc48-e65d-46dd-944c-ac41287da7dd
initial_sets, step_sizes =
	tm_integration(f, dom, k, J, time_horizon,
				   TIME_STEP_SIZE,
				   TIME_STEP_SIZE_EPS,
				   NR_CONTRACTIVENESS_TRIES,
				   NR_REFINEMENTS,
				   SCALE)

# ╔═╡ 9258bc2a-790d-4e3d-b250-b5a9ef0f0e81
md"""
# Plotting

## Preparation

Now we plot the results obtained from TM integration. Suppose we obtain the sequence of initial sets

$X_0, X_1, \dots, X_N$

then we first compute the sequence of corresponding boxes

$B_0, B_1, \dots, B_N$

These boxes constitute an over-approximation, an enclosure, of the initial sets.
"""

# ╔═╡ caf5cbb6-2613-4c1b-b13b-6affe34393f4
begin
	initial_boxes = map((Xi) -> interval_initial_set(Xi, dom), initial_sets)
	# FIXME: Assume t is the last variable, and drop its dummy vector element.
	initial_boxes = map((Bij) -> Bij[1:end-1], initial_boxes)
end

# ╔═╡ 1191638c-a1da-41df-81e7-66b7dd4b7191
begin
	include("scripts/clean_intervals_3_3_10.jl")  # Import most recent known boxes

	USE_LOCAL_HORIZON = false 	# Boolean flag, 
	known_boxes = copy(get_clean_boxes()) 	# Check intersection between computed and known boxes

	# Compare as many known VS computed boxes as possible, so that
	# you can always generate some plot even if the lengths mismatch.
	known_boxes = known_boxes[1:min(length(known_boxes), length(initial_boxes))]
end

# ╔═╡ 22867a28-ebb4-4d04-bf0e-2847ff774091
md"""
Next we define some configuration variables for the plotting section.
"""

# ╔═╡ 69977cd7-f7c5-4316-afbc-97c68f35f549
md"""
## Forward Euler

To visually validate the correctness of the computed and known results, we also **approximate the solution of the system of ODEs using forward Euler**. We can then plot those approximations for a subset of initial conditions $x_0 \in X_0$, to manually check that they are all enclosed by the boxes $B_0, B_1, \dots$
"""

# ╔═╡ 6b6dcfd2-6fe8-459e-8de7-661d46c6a20f
begin
	# Compute the forward Euler approximation, which yields
	# 	tseries = all time steps where an approximation was computed
	# 	vseries = [
	# 		the x variable approximations,
	# 	]
	euler_step_size = 0.001
	datax_max = time_horizon
	init_conditions = [0.1]
    tseries, vseries = euler(ode3310!, datax_max, euler_step_size, init_conditions)

	# Display forward Euler results legibly
	plots = []
	for var_idx in range(1, length(init_conditions))
		data_var = map((y) -> y[var_idx], vseries)
		plt = plot(tseries, data_var, xlabel="t", ylabel=vars_no_t[var_idx], legend=false)
		push!(plots, plt)
	end
	plot(plots..., layout=(1, length(plots)))
end

# ╔═╡ 7e91a680-0595-4dea-83c9-797f627ff014
md"""
## TM integration plots

And now generate some plots.

the "composed" figure plots each variable individually against time. For our ODE, that means **x against t**.

Also plot the approximate solution of the ODE obtained by **forward Euler** as a red curve.

"""

# ╔═╡ 2a5f17d6-3f91-4469-b6fb-dd469d57df65
begin
	# Compute the intersection of known and computed boxes
	intersect_boxes = map(((init, known),) -> intersect(IntervalBox(init), IntervalBox(known)), zip(initial_boxes, known_boxes))

	# For each pair (Bi, Gi) of computed box Bi and known box Gi, we want the intersection to be non-empty.
	intersect_not_empty = map((int) -> emptyinterval() != int, intersect_boxes)
    @assert all(intersect_not_empty)
end

# ╔═╡ ea0f2220-fa08-4160-9eb4-bfc9f912cba8
begin
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

    plt_intersect = plot_boxes_ND(intersect_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                         title="Bi ∩ Gi")

	# Collect all plots, which we will use as subplots
    figures = [plt_computed, plt_known, plt_both, plt_intersect]

	# Find the single ylim to enclose all plots.
    ylims_val = hull(Interval.(ylims.(figures))...)
    ylims_val = (ylims_val.lo, ylims_val.hi)
    xlims_val = hull(Interval.(xlims.(figures))...)

    plt_composed = plot(figures..., layout=(1, length(figures)),
                        ylims=ylims_val, xticks=[xlims_val.lo, mid(xlims_val), xlims_val.hi])

	# Add forward Euler curve
    datax = USE_LOCAL_HORIZON ? [0.0, map((_) -> euler_step_size, tseries[2:end])...] : tseries
    datay = map((vec) -> vec[1], vseries)

    plot_recursively!(plt_composed, datax, [datay], vars_no_t)

	plt_composed # Pluto cell output
end

# ╔═╡ Cell order:
# ╟─6b2d0d14-a438-11ef-11d9-e9ba654fb47e
# ╟─45cae1dc-ad20-4a5c-94ec-188b2e0380f2
# ╠═abf81117-2365-4a81-8116-73baf6f8caaf
# ╠═2580ed6a-98dc-427b-a5d8-344c5a9ed4ad
# ╟─774572a0-cfbe-4d4d-926e-c0964715608e
# ╠═e9eaca80-5d03-45d0-ace0-9e9939200aa6
# ╟─73c9ed6c-391c-4628-809b-60ec755afb94
# ╠═7ec7bc48-e65d-46dd-944c-ac41287da7dd
# ╟─9258bc2a-790d-4e3d-b250-b5a9ef0f0e81
# ╠═caf5cbb6-2613-4c1b-b13b-6affe34393f4
# ╟─22867a28-ebb4-4d04-bf0e-2847ff774091
# ╠═1191638c-a1da-41df-81e7-66b7dd4b7191
# ╟─69977cd7-f7c5-4316-afbc-97c68f35f549
# ╠═6b6dcfd2-6fe8-459e-8de7-661d46c6a20f
# ╟─7e91a680-0595-4dea-83c9-797f627ff014
# ╠═2a5f17d6-3f91-4469-b6fb-dd469d57df65
# ╠═ea0f2220-fa08-4160-9eb4-bfc9f912cba8

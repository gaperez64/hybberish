### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ bda06530-a1b0-11ef-04f3-01111228f3f7
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ 3a322105-09d3-4ae4-bca5-ba95bee26c2d
begin
	using Plots
	using TaylorSeries
	using IntervalArithmetic
	include("tm_integration.jl")
end

# ╔═╡ 0edfc1f9-9744-468e-bf3c-1b1f41b588b0
md"""
# TM integration

This notebook demonstrates the TM integration algorithm with hardcoded inputs.
You can play around with these inputs in the cells below, 

It is based on **Flowstar examples 3.3.2 and 3.3.6**.
"""

# ╔═╡ 520fa1bf-8e53-4edd-9d7d-03e54b9dfa90
md"""
## Environment

We first set up the environment, by including the necessary packages and files, etc.
"""

# ╔═╡ e00ceace-7b6d-4d4f-a489-2fdb0273eb1c
md"""
## Configure TM integration

We now set up the hardcoded inputs. The comments in the markdown cells and the in-line code comments will describe what input we use, and why.
"""

# ╔═╡ e3aa973e-2a96-45c7-bcb4-40758cc687be
# Full TM integration, based on Example 3.3.2 and Example 3.3.6
begin
	k = 3   # The TM arithmetic and truncation order
	NR_CONTRACTIVENESS_TRIES = 5
	NR_REFINEMENTS           = 1
	SCALE                    = 2.0
	TIME_STEP_SIZE = 0.02
	TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
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
end

# ╔═╡ b512424d-cd4a-41ec-94b7-28f9b74d21c0
md"""
We use the 2D system of ODEs specified in the flowstar paper examples 3.3.2 and 3.3.6:

x' = $(f[1])\
y' = $(f[2])

Several other parameters are in use.

* All TM arithmetic will be order k=$k.
* The finite time horizon is $\Delta$ = time\_horizon = $time_horizon.
* The fixed time step-size is $\delta$ = TIME\_STEP\_SIZE = $TIME_STEP_SIZE.
* Any step size $\delta_i$ < $TIME_STEP_SIZE_EPS = TIME\_STEP\_SIZE\_EPS simply gets skipped.
* Depending on the value of $\Delta$, the final time step $\delta_N$ may be smaller: $0 < \delta_N \leq \delta$
* J = $J is the initial remainder estimate, which gets refined during remainder refinement.
* Multiple parameters control the **remainder refinement** step of TM integration.
  * NR\_CONTRACTIVENESS\_TRIES is how many times we attempt widening $J_0$ to find a contractive remainder, before we give up.
  * NR\_REFINEMENTS is the number of refinements we apply to the contractive remainder $J_1 \subseteq J_0$.
  * SCALE is the scaling factor we apply when the contractiveness test fails. So $J_1 \nsubseteq J_0 \Rightarrow J_0 \leftarrow J_0 \cdot SCALE$
* * *dom* = $dom are the interval domains of the ODE variables x, y, t.
"""

# ╔═╡ 865b3373-1c90-4c52-9727-02a608327b29
md"""
## Run TM integration

iteration |   -   |      1     |   $\dots$  | N
    :---  | :---: |    :---:   |     :---:  | :---:
init set  | $X_0$ |    $X_1$   |   $\dots$  | $X_N$
step size |   0   | $\delta_1$ |   $\dots$  | $\delta_N$


Given the chosen configuration and inputs, run TM integration to compute
* ⠀$X_0, X_1, \dots, X_N$, the sequence of initial sets.
* ⠀$0, \delta_1, \dots, \delta_N$, the corresponding sequence of time steps. The first element is $0$ because it corresponds to $X_0$ which is the known, initial conditions at time $t = 0$.
"""

# ╔═╡ 123ee33c-ab75-4238-b9fd-2eed795f5c8b
initial_sets, step_sizes =
	tm_integration(f, dom, k, J, time_horizon,
				   TIME_STEP_SIZE,
				   TIME_STEP_SIZE_EPS,
				   NR_CONTRACTIVENESS_TRIES,
				   NR_REFINEMENTS,
				   SCALE)

# ╔═╡ b610fdbe-d21f-4fc0-bbc9-67c55bf385e6
md"""
# Plotting

## Preparation

Now we plot the results obtained from TM integration. Suppose we obtain the sequence of initial sets

$X_0, X_1, \dots, X_N$

then we first compute the sequence of corresponding boxes

$B_0, B_1, \dots, B_N$

These boxes constitute an over-approximation, an enclosure, of the initial sets.

Given an $m$-dimensional system of ODEs with variables $x_1, \dots, x_m$, initial set $X_i = (p(x_1, \dots, x_m), I)$ which is a Taylor model and $dom = (I_1, \dots, I_m)$ where $I_j$ is the interval domain of variable $x_j$.
Then you can compute box $B_j = p(I_1, \dots, I_m) + I$.
"""

# ╔═╡ 14934eac-2bea-4d70-afa5-43448e3281f6
# For j = 1..N compute the boxes Bj from the initial sets Xj
begin
	initial_boxes = map((Xi) -> interval_initial_set(Xi, dom), initial_sets)
	# FIXME: Assume t is the last variable, and drop its dummy vector element.
	initial_boxes = map((Bij) -> Bij[1:end-1], initial_boxes)
end

# ╔═╡ 244cb2cc-323b-43b9-a7fd-4ca2a1bff0d6
md"""
Next we define some configuration variables for the plotting section.

USE\_LOCAL\_HORIZON is a boolean flag. If **true**, then the $i$-th initial set will be plotted against the $i$-th **local** time step interval:

$[0, \delta_i]$

If **false**, then it is instead plotted against the $i$-th true (**global** so to say) time step interval

$[\delta_1 + \dots + \delta_{i-1}, \delta_1 + \dots + \delta_{i-1} + \delta_i] = \left[\sum_{j=1}^{i-1} \delta_j, \left( \sum_{j=1}^{i-1} \delta_j \right) + \delta_i\right]$
"""

# ╔═╡ 4e250e61-a17f-4742-a5fb-b5063562b421
begin
	USE_LOCAL_HORIZON = false 	# Boolean flag, 
	known_boxes = CLEAN_BOXES 	# Check intersection between computed and known boxes

	# Compare as many known VS computed boxes as possible, so that
	# you can always generate some plot even if the lengths mismatch.
	known_boxes = known_boxes[1:min(length(known_boxes), length(initial_boxes))]
end

# ╔═╡ a31e4811-a5cd-4656-a46a-eb8123522455
md"""
## Forward Euler

To visually validate the correctness of the computed and known results, we also **approximate the solution of the system of ODEs using forward Euler**. We can then plot those approximations for a subset of initial conditions $x_0 \in X_0$, to manually check that they are all enclosed by the boxes $B_0, B_1, \dots$
"""

# ╔═╡ 49709216-ee48-45f2-986e-3a191746d9bf
begin
	# Compute the forward Euler approximation, which yields
	# 	tseries = all time steps where an approximation was computed
	# 	vseries = [
	# 		the x variable approximations,
	# 		the y variable approximations
	# 	]
	euler_step_size = 0.001
	datax_max = time_horizon
	init_conditions = [0.0; 0.0]
	tseries, vseries = euler(ode!, datax_max, euler_step_size, init_conditions)

	# Display forward Euler results legibly
	plots = []
	for var_idx in range(1, length(init_conditions))
		data_var = map((y) -> y[var_idx], vseries)
		plt = plot(tseries, data_var, xlabel="t", ylabel=vars_no_t[var_idx], legend=false)
		push!(plots, plt)
	end
	plot(plots..., layout=(1, length(plots)))
end

# ╔═╡ 3dbb7e41-4fab-4b3e-8867-5d42c22765ad
md"""
## TM integration plots

And now generate some plots.

First, the "composed" figure plots each variable individually against time. For out ODE, that means **x against t, and y against t**. So, each of the rows corresponds to one ODE variable.

The plot is also divided into four columns:
* **Computed boxes Bi** : the boxes we derived above from the initial sets our TM integration function computed
* **Given boxes Gi** : the known boxes to compare to the computed boxes, e.g. the boxes computed by Flowstar
* **Bi & Gi** : combine the previous two columns in one plot, as a sanity check for the next column
* **Bi $\cap$ Gi** : the intersection of the computed and known boxes. i.e. for each pair Bi, Gi we compute their intersection Ii = Bi $\cap$ Gi and then draw Ii for all i

"""

# ╔═╡ 0ef5a8ed-051a-4faf-8768-68ddb19e8b93
begin
	# Compute the intersection of known and computed boxes
	intersect_boxes = map(((init, known),) -> intersect(IntervalBox(init), IntervalBox(known)), zip(initial_boxes, known_boxes))

	# For each pair (Bi, Gi) of computed box Bi and known box Gi, we want the intersection to be non-empty.
	intersect_not_empty = map((int) -> emptyinterval() != int, intersect_boxes)
    @assert all(intersect_not_empty)

	intersect_not_empty 	# Pluto cell output
end

# ╔═╡ 5922d398-8e55-4a0f-8a34-a44e7e9b024d
begin
	# Plot the sequence of computed boxes
    plt_computed = plot_boxes_ND(initial_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                         title="Computed boxes Bi", titlefontsize=8)

    # Plot the sequence of known boxes
    plt_known = plot_boxes_ND(known_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                         title="Given boxes Gi", titlefontsize=8)

	# Plot both sequences in one figure
	# Use a hack to plot TWO sequences of boxes on the same figure:
	# 	Append the box vectors to each other. BUT, a shifted version of the
	# 	time step vector, since the global time in the plotting function
	# 	increases for each box plotted.
    plt_both = nothing
    if USE_LOCAL_HORIZON
        plt_both = plot_boxes_ND([initial_boxes..., known_boxes...], [step_sizes..., step_sizes...], vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                            title="Bi & Gi")
    else
        plt_both = plot_boxes_ND([initial_boxes..., reverse(known_boxes)...], [step_sizes..., (-step_sizes)...], vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                            title="Bi & Gi")
    end

	# Plot only the intersection
    plt_intersect = plot_boxes_ND(intersect_boxes, step_sizes, vars_no_t, use_local_horizon=USE_LOCAL_HORIZON,
                         title="Bi ∩ Gi")

	# Collect all plots, which we will use as subplots
    figures = [plt_computed, plt_known, plt_both, plt_intersect]

    # Find the single ylim to enclose all plots' curves.
    ylims_val = hull(Interval.(ylims.(figures))...)
    ylims_val = (ylims_val.lo, ylims_val.hi)
    xlims_val = hull(Interval.(xlims.(figures))...)

    plt_composed = plot(figures..., layout=(1, length(figures)),
                        ylims=ylims_val, xticks=[xlims_val.lo, mid(xlims_val), xlims_val.hi])
end

# ╔═╡ 2f61baed-3aa8-42c3-846b-64512c95a8f3
md"""
Add forwar Euler to the composed plot, as a red line that describes the evolution of the plotted variable over time.
"""

# ╔═╡ e0147f92-6988-4d7a-b767-579a52a42961
begin
    datax_1 = USE_LOCAL_HORIZON ? [0.0, map((_) -> euler_step_size, tseries[2:end])...] : tseries
    datay_1 = map((vec) -> vec[2], vseries)

    plot_recursively!(plt_composed, datax_1, [datay_1, datay_1], vars_no_t)

	plt_composed # Pluto cell output
end

# ╔═╡ 8ec82d1b-eb8f-47d7-8764-c3a16fb26628
md"""
Alternatively, we could **plot the variables x and y against each other**, as opposed to plotting them individually against time t.
"""

# ╔═╡ ee8b5fb7-e9f5-4757-9a46-98e9a9772b45
begin
    plt_2D_computed =
		plot_boxes_2D(initial_boxes, vars_no_t, title="Computed boxes")
    plt_2D_known =
		plot_boxes_2D(known_boxes, vars_no_t, title="Known boxes")
    datax_2 = map((vec) -> vec[1], vseries)
    datay_2 = map((vec) -> vec[2], vseries)

    plot_recursively!(plt_2D_computed, datax_2, datay_2)
    plot_recursively!(plt_2D_known, datax_2, datay_2)

	# Compose the 2D plots
	plt_2D = plot(plt_2D_computed, plt_2D_known)
	plt_2D # Pluto cell output
end

# ╔═╡ Cell order:
# ╟─0edfc1f9-9744-468e-bf3c-1b1f41b588b0
# ╟─520fa1bf-8e53-4edd-9d7d-03e54b9dfa90
# ╠═bda06530-a1b0-11ef-04f3-01111228f3f7
# ╠═3a322105-09d3-4ae4-bca5-ba95bee26c2d
# ╟─e00ceace-7b6d-4d4f-a489-2fdb0273eb1c
# ╠═e3aa973e-2a96-45c7-bcb4-40758cc687be
# ╟─b512424d-cd4a-41ec-94b7-28f9b74d21c0
# ╟─865b3373-1c90-4c52-9727-02a608327b29
# ╠═123ee33c-ab75-4238-b9fd-2eed795f5c8b
# ╟─b610fdbe-d21f-4fc0-bbc9-67c55bf385e6
# ╠═14934eac-2bea-4d70-afa5-43448e3281f6
# ╟─244cb2cc-323b-43b9-a7fd-4ca2a1bff0d6
# ╠═4e250e61-a17f-4742-a5fb-b5063562b421
# ╟─a31e4811-a5cd-4656-a46a-eb8123522455
# ╠═49709216-ee48-45f2-986e-3a191746d9bf
# ╟─3dbb7e41-4fab-4b3e-8867-5d42c22765ad
# ╠═0ef5a8ed-051a-4faf-8768-68ddb19e8b93
# ╠═5922d398-8e55-4a0f-8a34-a44e7e9b024d
# ╟─2f61baed-3aa8-42c3-846b-64512c95a8f3
# ╠═e0147f92-6988-4d7a-b767-579a52a42961
# ╟─8ec82d1b-eb8f-47d7-8764-c3a16fb26628
# ╠═ee8b5fb7-e9f5-4757-9a46-98e9a9772b45

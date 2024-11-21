### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 3b649995-2db2-4e42-ad58-a5e95460d574
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ 33edc776-1ada-4fec-896a-426cd4614d74
begin
	using Plots
	using TaylorModels
	include("plotting.jl")
	include("euler.jl")
	include("tm_integration.jl")
	include("ode_shifting.jl")
end

# ╔═╡ 9ef48c3e-a529-11ef-3721-1330db209c3b
md"""
# TM integration

This notebook demonstrates the TM integration algorithm with hardcoded inputs.
You can play around with these inputs in the cells below, 

It is based on **Flowstar example 3.3.11**.
"""

# ╔═╡ f4f8b389-7090-4084-8d57-3e0e2907a27b
md"""
## Environment

We first set up the environment, by including the necessary packages and files, etc.
"""

# ╔═╡ 259509ae-a32b-4039-8bd0-d18a5dec163f
md"""
## ODE shift helper functions

We define some helper functions that aid in shifting/translating the ODEs, so that the initial hyperrectangle is centered on the origin.
"""

# ╔═╡ b879632b-4457-4f0f-b3d4-69ae5eb6e59f
md"""
## Configure TM integration

We now set up the hardcoded inputs. The comments in the markdown cells will be kept short. For more detailed explanations, see the notebook for example 3.3.6.
"""

# ╔═╡ f9d36075-e6bd-49bf-af79-0ce17108f282
begin
    # Full TM integration, based on Example 3.3.11
    k = 10   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 5
    NR_REFINEMENTS           = 10
    SCALE                    = 2.0
    TIME_STEP_SIZE = 0.02
    TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
    time_horizon = 3.1     # The finite time horizon
    vars = set_variables("x y t", order=k)
    vars_no_t = get_variable_names()[1:end-1]

    dom::IntervalBox = IntervalBox([
		4.9..5.1,   # x
		1.9..2.1,   # y
    	0..0.02  # t
	])

	# Shift the domains, and replace all variables by shifted variables.
	dom, vars_tms, shift_offsets = shifted_vars_tms(dom)

    # FIXME: Bumping up the order to avoid assertion errors
    # during TM arithmetic.
    old_order = get_order()
    new_order = old_order*2
    vars = set_variables(get_variable_string(), order=new_order)

    # The vector field f of the ODEs:
    #   f[1] = 1.5x - xy
    #   f[2] = -3y + xy
    # Represent the vector field as TMs:
    #   f[1] = 1.5x - xy + [0, 0]
    #   f[2] = -3y + xy + [0, 0]
    # where the Lagrange remainder is always zero, since
    # the vector field contains only polynomial terms
    f = [
		1.5 * vars_tms[1] - vars_tms[1] * vars_tms[2],
		-3 * vars_tms[2]  + vars_tms[1] * vars_tms[2]
	]
    # Initial remainder estimate J, a hyperrectangle
    J = fill(-0.1..0.1, length(f))

    # TODO: Reset to old order.
    vars = set_variables(get_variable_string(), order=old_order)

	"k=$k", "NR_REFINEMENTS=$NR_REFINEMENTS", "TIME_STEP_SIZE=$TIME_STEP_SIZE", "time_horizon=$time_horizon", dom, vars_tms, f # Pluto cell output
end

# ╔═╡ 45ede61b-c358-44c0-a096-9315fd0524fa
md"""
## Run TM integration

iteration |   -   |      1     |   $\dots$  | N
    :---  | :---: |    :---:   |     :---:  | :---:
init set  | $X_0$ |    $X_1$   |   $\dots$  | $X_N$
step size |   0   | $\delta_1$ |   $\dots$  | $\delta_N$
"""

# ╔═╡ d740e214-2303-4fc9-8cd6-bed4db6f98ad
initial_sets, step_sizes =
	tm_integration(f, dom, k, J, time_horizon,
				   TIME_STEP_SIZE,
				   TIME_STEP_SIZE_EPS,
				   NR_CONTRACTIVENESS_TRIES,
				   NR_REFINEMENTS,
				   SCALE)

# ╔═╡ b0b87d3d-c250-4eba-b3fe-6b6b193a9820
md"""
# Plotting

## Preparation

Now we plot the results obtained from TM integration. Suppose we obtain the sequence of initial sets

$X_0, X_1, \dots, X_N$

then we first compute the sequence of corresponding boxes

$B_0, B_1, \dots, B_N$

These boxes constitute an over-approximation, an enclosure, of the initial sets.
"""

# ╔═╡ 1d6e06a1-43ad-4ba7-b89e-598e72620b35
begin
	initial_boxes = map((Xi) -> interval_initial_set(Xi, dom), initial_sets)
	# FIXME: Assume t is the last variable, and drop its dummy vector element.
	initial_boxes = map((Bij) -> Bij[1:end-1], initial_boxes)
	# Undo the shift applied to the ODEs and domains, to obtain the
	# boxes of the true ODEs and domains.
	initial_boxes = unshift_boxes(initial_boxes, shift_offsets)
	# Typecast the initial boxes for cleanliness of the source code.
	initial_boxes = map((Bij) -> IntervalBox(Bij...), initial_boxes)
end

# ╔═╡ f048f8ad-feea-4603-9ecc-889eb5ce288a
begin
	include("scripts/clean_intervals_3_3_11.jl")  # Import most recent known boxes

	# Boolean flags to toggle which subplots to show and how.
	SHOW_ONLY_COMPUTED_BOXES = true
	USE_LOCAL_HORIZON = false

	# The number of flowpipes/boxes to keep.
	trunc_max = 40
	trunc_max = length(initial_boxes)

	# Check intersection between computed and known boxes
	known_boxes = copy(get_clean_boxes())
	known_boxes = map((_) -> known_boxes[1], eachindex(initial_boxes))

	# Compare as many known VS computed boxes as possible, so that
	# you can always generate some plot even if the lengths mismatch.
	known_boxes = known_boxes[1:min(length(known_boxes), length(initial_boxes))]

	# FIXME: just ignore known boxes since it's dummy
	empty!(known_boxes)

	length_diff = length(initial_boxes) - length(known_boxes)
	if length_diff > 0
		padding = map((_) -> IntervalBox(zero(initial_boxes[1])), 1:length_diff)
		known_boxes = vcat(known_boxes, padding)
	end
end

# ╔═╡ 84afb4c5-c7b4-40fd-bb1a-37ab44994eb9
md"""
Next we define some configuration variables for the plotting section.
"""

# ╔═╡ 3ede02f6-2b5f-4476-8c6d-f1f4a9b275ca
md"""
## TM integration plots

And now generate some plots.

the "composed" figure plots each variable individually against time. For our ODE, that means **x against t**.

Also plot the approximate solution of the ODE obtained by **forward Euler** as a red curve.

"""

# ╔═╡ 960db067-b6bd-41a1-a1d9-be30d42ad2ae
begin
	# Truncate several series, so that you can play around with the number
	# of flowpipes/boxes to display on the plots without having to re-run
	# TM integration.

	# Truncating should be safe, only the plots follow after this.
	time_horizon_truncated = trunc_max * TIME_STEP_SIZE
	step_sizes_truncated = step_sizes[1:trunc_max]
	initial_boxes_truncated = initial_boxes[1:trunc_max]
	known_boxes_truncated = known_boxes[1:trunc_max]
end

# ╔═╡ 87acf7ba-fd57-40e0-b9e7-53ff97379a1f
begin
	euler_init_state = [
		5.0; # x
		2.0  # y
	]
	# Generate the composed plot, which also does euler
	plt_composed, tseries_, vseries_ = plot_vars_against_time(initial_boxes_truncated, known_boxes_truncated, time_horizon_truncated, step_sizes_truncated, vars_no_t, USE_LOCAL_HORIZON, ode3311!, euler_init_state)

	if SHOW_ONLY_COMPUTED_BOXES
		plot(plot(plt_composed[1]), plot(plt_composed[2]), layout=(2,1)) # Pluto cell output
	else
		plt_composed # Pluto cell output
	end
end

# ╔═╡ 3438d00a-7ccd-4599-83d0-1eec4c3330db
begin
	plt_2D_computed =
		plot_boxes_2D(initial_boxes_truncated, vars_no_t, title="Computed boxes")
    plt_2D_known =
		plot_boxes_2D(known_boxes_truncated, vars_no_t, title="Known boxes")
    datay_2x = map((vec) -> vec[1], vseries_)
    datay_2y = map((vec) -> vec[2], vseries_)

    plot_recursively!(plt_2D_computed, datay_2x, datay_2y)
    plot_recursively!(plt_2D_known, datay_2x, datay_2y)

	# Compose the 2D plots
	plt_2D = plot(plt_2D_computed, plt_2D_known)

	if SHOW_ONLY_COMPUTED_BOXES
		plot(plt_2D[1]) # Pluto cell output
	else
		plt_2D # Pluto cell output
	end
end

# ╔═╡ Cell order:
# ╟─9ef48c3e-a529-11ef-3721-1330db209c3b
# ╟─f4f8b389-7090-4084-8d57-3e0e2907a27b
# ╠═3b649995-2db2-4e42-ad58-a5e95460d574
# ╠═33edc776-1ada-4fec-896a-426cd4614d74
# ╟─259509ae-a32b-4039-8bd0-d18a5dec163f
# ╟─b879632b-4457-4f0f-b3d4-69ae5eb6e59f
# ╠═f9d36075-e6bd-49bf-af79-0ce17108f282
# ╟─45ede61b-c358-44c0-a096-9315fd0524fa
# ╠═d740e214-2303-4fc9-8cd6-bed4db6f98ad
# ╟─b0b87d3d-c250-4eba-b3fe-6b6b193a9820
# ╠═1d6e06a1-43ad-4ba7-b89e-598e72620b35
# ╟─84afb4c5-c7b4-40fd-bb1a-37ab44994eb9
# ╠═f048f8ad-feea-4603-9ecc-889eb5ce288a
# ╟─3ede02f6-2b5f-4476-8c6d-f1f4a9b275ca
# ╠═960db067-b6bd-41a1-a1d9-be30d42ad2ae
# ╠═87acf7ba-fd57-40e0-b9e7-53ff97379a1f
# ╠═3438d00a-7ccd-4599-83d0-1eec4c3330db

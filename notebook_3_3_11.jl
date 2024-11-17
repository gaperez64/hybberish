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

# ╔═╡ b879632b-4457-4f0f-b3d4-69ae5eb6e59f
md"""
## Configure TM integration

We now set up the hardcoded inputs. The comments in the markdown cells will be kept short. For more detailed explanations, see the notebook for example 3.3.6.
"""

# ╔═╡ f9d36075-e6bd-49bf-af79-0ce17108f282
begin
    # Full TM integration, based on Example 3.3.11
    k = 3   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 5
    NR_REFINEMENTS           = 1
    SCALE                    = 2.0
    TIME_STEP_SIZE = 0.02
    TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
    time_horizon = 0.7     # The finite time horizon
    vars = set_variables("x y t", order=k)
    vars_no_t = get_variable_names()[1:end-1]

	# FIXME: offset x by +5 and y by +2, since the TaylorModels
	# library requires 0 to be in the domains and such.
    dom = IntervalBox([
		-1..1,   # x
		-1..1,   # y
    	0..0.02 # t
	])
    vars_tms = tm_initial_set(dom)
    
    
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
	# FIXME: offset x by +5 and y by +2, since the TaylorModels
	# library requires 0 to be in the domains and such.
    f = [
		1.5 * vars_tms[1] - vars_tms[1] * vars_tms[2] + 5,
		-3 * vars_tms[2]  + vars_tms[1] * vars_tms[2] + 2
	]
    # Initial remainder estimate J, a hyperrectangle
    J = fill(-0.1..0.1, length(f))

    # TODO: Reset to old order.
    vars = set_variables(get_variable_string(), order=old_order)

	f # Pluto cell output
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
	# Typecast the initial boxes for cleanliness of the source code.
	initial_boxes = map((Bij) -> IntervalBox(Bij...), initial_boxes)
end

# ╔═╡ f048f8ad-feea-4603-9ecc-889eb5ce288a
begin
	include("scripts/clean_intervals_3_3_11.jl")  # Import most recent known boxes

	USE_LOCAL_HORIZON = false 	# Boolean flag, 
	known_boxes = copy(get_clean_boxes()) 	# Check intersection between computed and known boxes

	# Compare as many known VS computed boxes as possible, so that
	# you can always generate some plot even if the lengths mismatch.
	known_boxes = known_boxes[1:min(length(known_boxes), length(initial_boxes))]
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

# ╔═╡ 87acf7ba-fd57-40e0-b9e7-53ff97379a1f
begin
	euler_init_state = [
		0.0; # x
		0.0  # y
	]
	# Generate the composed plot, which also does euler
	plt_composed, tseries_, vseries_ = plot_vars_against_time(initial_boxes, known_boxes, time_horizon, step_sizes, vars_no_t, USE_LOCAL_HORIZON, ode332!, euler_init_state)

	plt_composed # Pluto cell output
end

# ╔═╡ Cell order:
# ╟─9ef48c3e-a529-11ef-3721-1330db209c3b
# ╟─f4f8b389-7090-4084-8d57-3e0e2907a27b
# ╠═3b649995-2db2-4e42-ad58-a5e95460d574
# ╠═33edc776-1ada-4fec-896a-426cd4614d74
# ╟─b879632b-4457-4f0f-b3d4-69ae5eb6e59f
# ╠═f9d36075-e6bd-49bf-af79-0ce17108f282
# ╟─45ede61b-c358-44c0-a096-9315fd0524fa
# ╠═d740e214-2303-4fc9-8cd6-bed4db6f98ad
# ╟─b0b87d3d-c250-4eba-b3fe-6b6b193a9820
# ╠═1d6e06a1-43ad-4ba7-b89e-598e72620b35
# ╟─84afb4c5-c7b4-40fd-bb1a-37ab44994eb9
# ╠═f048f8ad-feea-4603-9ecc-889eb5ce288a
# ╟─3ede02f6-2b5f-4476-8c6d-f1f4a9b275ca
# ╠═87acf7ba-fd57-40e0-b9e7-53ff97379a1f

### A Pluto.jl notebook ###
# v0.20.4

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
	include("../src/old/tm_integration.jl")
	include("../src/euler.jl")
end

# ╔═╡ 039988bf-1261-495d-b3de-be812378168f
include("../src/plotting.jl")

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

# ╔═╡ b1869321-9073-4ec7-b4ea-d30ae37bb43e
"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function vector_field!(du::Vector, u::Vector)
	x, y = u

	du[1] = 1 + y
	du[2] = -(x^2)
end

# ╔═╡ e3aa973e-2a96-45c7-bcb4-40758cc687be
# Full TM integration, based on Example 3.3.2 and Example 3.3.6
begin
	k = 10   # The TM arithmetic and truncation order
	NR_CONTRACTIVENESS_TRIES = 10
	NR_REFINEMENTS           = 10
	SCALE                    = 2.0
	TIME_STEP_SIZE = 0.01
	TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
	time_horizon = 0.02     # The finite time horizon
	vars = set_variables("x y t", order=k)

	dom = IntervalBox([
		-1..1,      # x
		-0.5..0.5,  # y
	])

	# Initial remainder estimate J, a hyperrectangle
    J = IntervalBox(fill(-0.1..0.1, length(dom)))

	PRECONDITIONING::Bool = true
	RECOMPUTE_VECTOR_FIELD::Bool = false
end

# ╔═╡ b512424d-cd4a-41ec-94b7-28f9b74d21c0
md"""
We use the 2D system of ODEs specified in the flowstar paper examples 3.3.2 and 3.3.6.

Several parameters are in use.

* All TM arithmetic will be order k=$k.
* The finite time horizon is $\Delta$ = time\_horizon = $time_horizon.
* The fixed time step-size is $\delta$ = TIME\_STEP\_SIZE = $TIME_STEP_SIZE.
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
begin
	over_approximations::Vector{OverApprox} =
		tm_integration(vector_field!, dom, k, J, time_horizon,
					   TIME_STEP_SIZE,
					   TIME_STEP_SIZE_EPS,
					   NR_CONTRACTIVENESS_TRIES,
					   NR_REFINEMENTS,
					   SCALE,
					   PRECONDITIONING=PRECONDITIONING,
		               RECOMPUTE_VECTOR_FIELD=RECOMPUTE_VECTOR_FIELD)

	initial_sets = map((approx) -> approx.initial_set, over_approximations)
	step_sizes = map((approx) -> approx.step_size, over_approximations)
	dynamic_boxes = map((approx) -> approx.updated_box, over_approximations)
	initial_boxes = map((approx) -> approx.initial_box, over_approximations)

	initial_sets, step_sizes, initial_boxes, dynamic_boxes # Pluto cell output
end

# ╔═╡ b610fdbe-d21f-4fc0-bbc9-67c55bf385e6
md"""
# Plotting

## Preparation

Next we define some configuration variables for the plotting section.

USE\_LOCAL\_HORIZON is a boolean flag. If **true**, then the $i$-th initial set will be plotted against the $i$-th **local** time step interval:

$[0, \delta_i]$

If **false**, then it is instead plotted against the $i$-th true (**global** so to say) time step interval

$[\delta_1 + \dots + \delta_{i-1}, \delta_1 + \dots + \delta_{i-1} + \delta_i] = \left[\sum_{j=1}^{i-1} \delta_j, \left( \sum_{j=1}^{i-1} \delta_j \right) + \delta_i\right]$
"""

# ╔═╡ 3c240a9f-2b22-4b88-b4d1-50b96056adad
"""A representation of the following 2D ODE:
    x' = 1 + y
    y' = -x^2
"""
function euler!(du, u, p, t)
    du[1] = 1.0 + u[2]
    du[2] = -(u[1]^2)
end

# ╔═╡ ce9d4671-f25f-4176-8818-7f0a193643a8
begin
    vars_no_t = get_variable_names()[1:end-1]
	euler_init_state = [
		0.0,  # x
		0.0,  # x
	]

    tseries, vseries = euler(euler!, time_horizon, TIME_STEP_SIZE, euler_init_state)
end

# ╔═╡ 4f43995b-8286-49f1-b2f1-2b666961100b
begin
	PLOT_DYNAMIC_BOXES::Bool = false

	if PLOT_DYNAMIC_BOXES
		to_draw_boxes = dynamic_boxes
	else
		to_draw_boxes = initial_boxes
	end
end

# ╔═╡ e13d8852-828d-455e-830d-a6bd60640782
begin
	pltND = plot_boxes_ND(to_draw_boxes, step_sizes, vars_no_t)

	# for 

	pltND # Pluto cell output
end

# ╔═╡ 24c38d62-ff4d-45f7-b671-0705c7ac752f
begin
	# Generate the composed plot, which also does euler
	plt2D = plot_boxes_2D(to_draw_boxes, Tuple(vars_no_t))

	plt2D # Pluto cell output
end

# ╔═╡ Cell order:
# ╟─0edfc1f9-9744-468e-bf3c-1b1f41b588b0
# ╟─520fa1bf-8e53-4edd-9d7d-03e54b9dfa90
# ╠═bda06530-a1b0-11ef-04f3-01111228f3f7
# ╠═3a322105-09d3-4ae4-bca5-ba95bee26c2d
# ╟─e00ceace-7b6d-4d4f-a489-2fdb0273eb1c
# ╠═b1869321-9073-4ec7-b4ea-d30ae37bb43e
# ╠═e3aa973e-2a96-45c7-bcb4-40758cc687be
# ╟─b512424d-cd4a-41ec-94b7-28f9b74d21c0
# ╟─865b3373-1c90-4c52-9727-02a608327b29
# ╠═123ee33c-ab75-4238-b9fd-2eed795f5c8b
# ╟─b610fdbe-d21f-4fc0-bbc9-67c55bf385e6
# ╠═039988bf-1261-495d-b3de-be812378168f
# ╠═3c240a9f-2b22-4b88-b4d1-50b96056adad
# ╠═ce9d4671-f25f-4176-8818-7f0a193643a8
# ╠═4f43995b-8286-49f1-b2f1-2b666961100b
# ╠═e13d8852-828d-455e-830d-a6bd60640782
# ╠═24c38d62-ff4d-45f7-b671-0705c7ac752f

### A Pluto.jl notebook ###
# v0.20.1

using Markdown
using InteractiveUtils

# ╔═╡ fd505954-161d-4854-adb7-dfbcf5bdd076
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate(".")
end

# ╔═╡ 45a71bbe-6f74-4453-9b89-ec9621848b95
begin
	using Plots
	using TaylorModels
	include("plotting.jl")
	include("euler.jl")
	include("tm_integration.jl")
end

# ╔═╡ 28be9778-a531-11ef-225c-691fbd9f45d1
md"""
# TM integration

This notebook demonstrates the TM integration algorithm with hardcoded inputs.
You can play around with these inputs in the cells below, 

It is based on **the truck trailer example.**.
"""

# ╔═╡ 2db4175f-ef9a-42ba-a4f0-4b6bedce35df
md"""
## Environment

We first set up the environment, by including the necessary packages and files, etc.
"""

# ╔═╡ ed41723e-29e4-4f2b-90af-9565ff83f061
md"""
## Configure TM integration

We now set up the hardcoded inputs. The comments in the markdown cells will be kept short. For more detailed explanations, see the notebook for example 3.3.6.
"""

# ╔═╡ 5f7da750-ee2c-4f69-9a65-2314f10cd91b
"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function vector_field!(du::Vector, u::Vector)
	x1, y1, theta0, theta1, v0, dtheta0, = u

    L0 = 0.3375  # truck length
    M0 = 0.1     # truck distance to center of mass
    L1 = 0.3     # trailer length
    M1 = 0.06    # trailer distance to center of mass

	# intermediate variables
	beta01 = theta0 - theta1  # angle between truck and trailer
	v1 = v0 * cos(beta01) + M0 * sin(beta01) * dtheta0  # trailer velocity

	# ODEs + zero(u[.]) for type stability
	du[1] = v1 * cos(theta1) + zero(u[1])
	du[2] = v1 * sin(theta1) + zero(u[2])
	du[3] = dtheta0 + zero(u[3])
	du[4] = (v0 / L1) * sin(beta01) - (M0 / L1) * cos(beta01) * dtheta0 + zero(u[4])
	du[5] = zero(u[5])  # v0 remains constant
	du[6] = zero(u[6])  # dtheta0 remains constant
end

# ╔═╡ 1cc16d32-a6e6-4e4b-b2b2-f913c8cbc814
begin
    k = 2   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 10
    NR_REFINEMENTS           = 10
    SCALE                    = 2.0
    TIME_STEP_SIZE = 0.02
    TIME_STEP_SIZE_EPS = 2.0e-8  # FIXME: REMOVE THIS FROM THE SOURCE CODE
    time_horizon = 0.2     # The finite time horizon
    vars = set_variables("x y θ0 θ1 v0 dθ0 t", order=k)

    dom = IntervalBox([
		0..0,   # x
		0.25..0.25,   # y
		1.57..1.57,   # θ0
		1.57..1.57,   # θ1
		0.1..0.3,   # v0
		0.2..0.4,   # dθ0
	])

	# Initial remainder estimate J, a hyperrectangle
	J = IntervalBox(fill(-0.1..0.1, length(dom)))

	PRECONDITIONING::Bool = true
	RECOMPUTE_VECTOR_FIELD::Bool = false

	"k=$k", "NR_REFINEMENTS=$NR_REFINEMENTS", "TIME_STEP_SIZE=$TIME_STEP_SIZE", "time_horizon=$time_horizon", "domains=$dom", "typeof(dom)=$(typeof(dom))" # Pluto cell output
end

# ╔═╡ 6f6fcac2-fdc4-41df-8f9e-2c2640ff2228
md"""
## Run TM integration

iteration |   -   |      1     |   $\dots$  | N
    :---  | :---: |    :---:   |     :---:  | :---:
init set  | $X_0$ |    $X_1$   |   $\dots$  | $X_N$
step size |   0   | $\delta_1$ |   $\dots$  | $\delta_N$
"""

# ╔═╡ a90b2b22-cb32-4a54-b45c-1e926d99e8ae
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

# ╔═╡ 130e6a05-936b-44e0-b456-8979f552d714
md"""
# Plotting

## Preparation

Now we plot the results obtained from TM integration.
"""

# ╔═╡ 59a2716b-ebf7-4cc8-9e43-f0817a613827
begin
	# If true, then plots the boxes computed using the domains that are recomputed every iteration of TM integration.
	# If false, then plots the boxes computed using exclusively the initial domains/hyperrectangle.
	PLOT_DYNAMIC_BOXES::Bool = false

	draw_boxes = PLOT_DYNAMIC_BOXES ? dynamic_boxes : initial_boxes
end

# ╔═╡ c3f65f57-afc4-49b2-9bd4-5e11f9c9f1f1
md"""
Next we define some configuration variables for the plotting section.
"""

# ╔═╡ 66ef3873-c5c3-4d0e-b6ea-58681ae49c1e
md"""
## TM integration plots

And now generate some plots.

We also plot the approximate solution of the ODE obtained by **forward Euler** as a red curve. So we first need to evaluate the ODEs with Forward Euler first.

"""

# ╔═╡ ecaba1a6-6a5d-4d86-85ea-f893b4557290
"""The Forward Euler vector field."""
function vector_field_euler!(du, u, p, t)
	return vector_field!(du, u)
end

# ╔═╡ 8c18009b-de7b-4ecf-94a3-f8f874035628
begin
    vars_no_t = get_variable_names()[1:end-1]

	# WARNING: this is overwritten below by default!!!!!!!!!!!!!
	euler_init_state = [
		0.0;  # x
		0.25; # y
		1.57; # θ0
		1.57; # θ1
		0.2;  # v0
		0.3   # dθ0
	]
	# Start Forward Euler in the middle of the variable domains.
	# TODO: Comment this line if you want to manually control u0.
	euler_init_state = Vector(mid(dom))

	# Evaluate Forward Euler.
    eseries = euler(vector_field_euler!, time_horizon, TIME_STEP_SIZE, euler_init_state)
end

# ╔═╡ fcfa1c8c-8807-4d80-8dfc-4a15002047e8
begin
	# Actual plotting
	pltND = plot_boxes_ND(draw_boxes, step_sizes, vars_no_t)

	plotsND = [plot(p) for p in pltND.subplots] # Pluto cell output
	plot!.(plotsND, eseries)
end

# ╔═╡ 385c7d3b-5c47-4bc9-b03b-b1fd3aa9b157
begin
	# Extract the x,y components from each box.
	xy_boxes = [
		IntervalBox(box[1:2])
		for box in draw_boxes
	]

	# Extract the approximated values for the variables x and y.
	# Zip the two lists, to obtain a list [(x, y), (x, y), ...].
	xy_forward_euler_data = collect(zip(
		getindex.(eseries[1], 2),
		getindex.(eseries[2], 2)
	))
	# Actual plotting.
	plt2D = plot_boxes_2D(xy_boxes, (vars_no_t[1], vars_no_t[2]))
	plot!(plt2D, xy_forward_euler_data)
end

# ╔═╡ Cell order:
# ╟─28be9778-a531-11ef-225c-691fbd9f45d1
# ╟─2db4175f-ef9a-42ba-a4f0-4b6bedce35df
# ╠═fd505954-161d-4854-adb7-dfbcf5bdd076
# ╠═45a71bbe-6f74-4453-9b89-ec9621848b95
# ╟─ed41723e-29e4-4f2b-90af-9565ff83f061
# ╠═5f7da750-ee2c-4f69-9a65-2314f10cd91b
# ╠═1cc16d32-a6e6-4e4b-b2b2-f913c8cbc814
# ╟─6f6fcac2-fdc4-41df-8f9e-2c2640ff2228
# ╠═a90b2b22-cb32-4a54-b45c-1e926d99e8ae
# ╟─130e6a05-936b-44e0-b456-8979f552d714
# ╠═59a2716b-ebf7-4cc8-9e43-f0817a613827
# ╟─c3f65f57-afc4-49b2-9bd4-5e11f9c9f1f1
# ╟─66ef3873-c5c3-4d0e-b6ea-58681ae49c1e
# ╠═ecaba1a6-6a5d-4d86-85ea-f893b4557290
# ╠═8c18009b-de7b-4ecf-94a3-f8f874035628
# ╠═fcfa1c8c-8807-4d80-8dfc-4a15002047e8
# ╠═385c7d3b-5c47-4bc9-b03b-b1fd3aa9b157

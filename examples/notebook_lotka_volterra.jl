### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 3cf832b2-0ef6-45db-91d6-d02cf0e767e9
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ bc1a46bf-71ca-461c-82f4-20b58dc9f3e2
begin
    include("../src/plotting.jl")
    include("../src/euler.jl")
    include("../src/tm_integration.jl")
	using IntervalArithmetic
	using TaylorSeries
	using Plots
end

# ╔═╡ 24f5b570-029e-11f0-2f76-ef0008d0b5f1
md"""
# Setup
"""

# ╔═╡ 215f1d90-d451-4de9-8fd4-6de2d4a9c018
md"""
# Define the Dynamics

We define a constructor function that is used to
1. Taylorize the dynamics as input to the Lie derivative routine.
2. Perform the TM composition in the TM extension of the picard operator.
"""

# ╔═╡ 99a28c51-3bd1-40ee-b59f-d90f25197145
"""Construct the dynamics.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function f_dot!(du::Vector, u::Vector)
    x, y, t = u
    du[1] = 1.5*x - x*y
	du[2] =  -3*y + x*y
end


# ╔═╡ 99cd7fd8-7af9-43d1-8f04-b79f4aea79d1
begin
	#
    # Full TM integration, based on Example 3.3.11 (lotka volterra systems)
	#

	# The truncation degree / the degree of all polynomials
	# that are used during computations.
	k = 10

	# The fixed time step size.
	tstep::Float64 = 0.01
	# Specify time as a finite time horizon.
	time_horizon::Float64 = 0.08
	
	# The scale factor for when contractiveness fails.
	nr_contractiveness_tries = 10
	# The number of refinements to perform at most.
	nr_refinements = 10
	# Quit refinement early if the improvement a single refinement
	# provides falls below this threshold.
	refinement_eps = 0.001
	# The scale factor with which to widen the initial safe remainder
	# estimate when the contractiveness check fails.
	scale = 2.0


	initial = [
		( "x", 4.9..5.1 ),
		( "y", 1.9..2.1 ),
		( "t", 0.0..0.0 )
    ]
    init_values = IntervalBox([e[2] for e in initial])

end

# ╔═╡ 16441a15-cf51-4d9e-aa80-20ada1e0b0ed
begin

	boxes::Vector{IntervalBox} = []
	fboxes::Vector{IntervalBox} = []
	boxes, fboxes = tm_integration(
	    f_dot!,
	    initial,
	    k,
	    (-0.1..0.1),
	    time_horizon,
	    tstep,
	    nr_contractiveness_tries,
	    nr_refinements,
	    SCALE=scale,
	    REFINEMENT_EPS=refinement_eps
	)
	vnames = get_variable_names()

end

# ╔═╡ 615e5677-445b-4dca-990d-cb61fa0a04b9
md"""
# Plotting

**If you want to truncate the series of boxes to inspect a smaller time horizon, see the cell that displays the plots below.**

"""

# ╔═╡ 603bb58e-4cb5-4351-8c67-582104c6dfa3
"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function ode_euler!(du, u, p, t)
	x, y = u
    du[1] = 1.5*x - x*y
	du[2] =  -3*y + x*y
end

# ╔═╡ f3bb1ee8-e06b-4897-a652-38a209e94a59
begin

	# Start Forward Euler in the middle of the variable domains.
	euler_init_state = Vector(mid(init_values))

	# Evaluate Forward Euler.
	euler_step_ratio::Integer = 10
	euler_step::Float64 = tstep / euler_step_ratio
	eseries = euler(ode_euler!, time_horizon, euler_step, euler_init_state)
end

# ╔═╡ 2a069428-d55f-4866-b7b3-5a119de5e795
begin

	# (Optional) Truncate the vector of boxes, to inspect
	# the TM integration results.
	nr_boxes = length(boxes) # TODO: Comment this line if you want to truncate
	nr_boxes = 8

	@assert nr_boxes <= length(boxes)
	@assert nr_boxes <= length(fboxes)

	truncated_boxes = boxes[1:nr_boxes]
	truncated_fboxes = fboxes[1:nr_boxes]
	truncated_eseries = [
		series[1:nr_boxes*euler_step_ratio]
		for series in eseries
	]

	tdata = [e[1] for e in truncated_eseries[1]]
	xdata = [e[2] for e in truncated_eseries[1]]
	ydata = [e[2] for e in truncated_eseries[2]]

	# Generate the plots
	pltND1 = plot_boxes_ND(truncated_boxes, vnames, sgtitle="vals", legend=true)
	plot!(pltND1.subplots[1], tdata, xdata, label="ODE Forward Euler")
	plot!(pltND1.subplots[2], tdata, ydata, label="ODE Forward Euler")

	pltND2 = plot_boxes_ND(truncated_fboxes, vnames, sgtitle="fpipe", legend=true)
	plot!(pltND2.subplots[1], tdata, xdata, label="ODE Forward Euler")
	plot!(pltND2.subplots[2], tdata, ydata, label="ODE Forward Euler")

	# Make Pluto notebook display the plots
	pltND = plot(pltND1, pltND2)

end

# ╔═╡ Cell order:
# ╟─24f5b570-029e-11f0-2f76-ef0008d0b5f1
# ╠═3cf832b2-0ef6-45db-91d6-d02cf0e767e9
# ╠═bc1a46bf-71ca-461c-82f4-20b58dc9f3e2
# ╟─215f1d90-d451-4de9-8fd4-6de2d4a9c018
# ╠═99a28c51-3bd1-40ee-b59f-d90f25197145
# ╠═99cd7fd8-7af9-43d1-8f04-b79f4aea79d1
# ╠═16441a15-cf51-4d9e-aa80-20ada1e0b0ed
# ╟─615e5677-445b-4dca-990d-cb61fa0a04b9
# ╠═603bb58e-4cb5-4351-8c67-582104c6dfa3
# ╠═f3bb1ee8-e06b-4897-a652-38a209e94a59
# ╠═2a069428-d55f-4866-b7b3-5a119de5e795

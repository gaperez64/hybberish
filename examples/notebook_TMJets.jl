### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ e0bbb60b-b654-42ca-8b18-f108f17669f4
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ 158e8e41-4f2d-4eff-aff1-fdee7c6d9515
begin
    include("../src/plotting.jl")
    include("../src/euler.jl")
    include("../src/tm_integration.jl")
	using IntervalArithmetic
	using TaylorSeries
	using Plots
end

# ╔═╡ 143689d0-02b4-11f0-28c9-770685f0b21d
md"""
# Setup
"""

# ╔═╡ e92b89a1-18bf-4c9e-b9af-d45bf1fe937d
md"""
# Define the Dynamics

We define a constructor function that is used to
1. Taylorize the dynamics as input to the Lie derivative routine.
2. Perform the TM composition in the TM extension of the picard operator.
"""

# ╔═╡ 6d61f50b-8cd7-4fb0-b06b-e4989d3cbe0b
"""Construct the dynamics.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function f_dot!(du::Vector, u::Vector)
	y, = u
	du[1] = y * cos(y)
end


# ╔═╡ 97a4186b-1989-4a20-b60e-b86bb6df9a6a
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
	time_horizon::Float64 = 0.10
	
	# The scale factor for when contractiveness fails.
	nr_contractiveness_tries = 100
	# The number of refinements to perform at most.
	nr_refinements = 100
	# Quit refinement early if the improvement a single refinement
	# provides falls below this threshold.
	refinement_eps = 0.001
	# The scale factor with which to widen the initial safe remainder
	# estimate when the contractiveness check fails.
	scale = 2.0


	initial = [
		( "x", 1.4..1.6 ),
		( "t", 0.0..0.0 )
    ]
    init_values = IntervalBox([e[2] for e in initial])

end

# ╔═╡ 76781e85-099c-4b2c-bb6e-96536523750e
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

# ╔═╡ 74f82d5c-7dbf-4aa6-b021-3a3881d9ad5c
md"""
# Plotting

**If you want to truncate the series of boxes to inspect a smaller time horizon, see the cell that displays the plots below.**

"""

# ╔═╡ 1963212f-6dbc-4ca1-aa8c-0f9958a5909d
"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function ode_euler!(du, u, p, t)
	y, = u
	du[1] = y * cos(y)
end

# ╔═╡ 74ccb6be-dd28-45f0-a731-f69a75a17fc3
begin

	# Start Forward Euler in the middle of the variable domains.
	euler_init_state = Vector(mid(init_values))

	# Evaluate Forward Euler.
	euler_step_ratio::Integer = 10
	euler_step::Float64 = tstep / euler_step_ratio
	eseries = euler(ode_euler!, time_horizon, euler_step, euler_init_state)
end

# ╔═╡ e11c2e4e-f27d-4813-9861-fa2a6f23e523
begin

	# (Optional) Truncate the vector of boxes, to inspect
	# the TM integration results.
	nr_boxes = 8
	nr_boxes = length(boxes) # TODO: Comment this line if you want to truncate

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

	# Generate the plots
	pltND1 = plot_boxes_ND(truncated_boxes, vnames, sgtitle="vals", legend=true)
	plot!(pltND1, tdata, xdata, label="ODE Forward Euler")

	pltND2 = plot_boxes_ND(truncated_fboxes, vnames, sgtitle="fpipe", legend=true)
	plot!(pltND2, tdata, xdata, label="ODE Forward Euler")

	# Make Pluto notebook display the plots
	pltND = plot(pltND1, pltND2)

end

# ╔═╡ Cell order:
# ╠═143689d0-02b4-11f0-28c9-770685f0b21d
# ╠═e0bbb60b-b654-42ca-8b18-f108f17669f4
# ╠═158e8e41-4f2d-4eff-aff1-fdee7c6d9515
# ╟─e92b89a1-18bf-4c9e-b9af-d45bf1fe937d
# ╠═6d61f50b-8cd7-4fb0-b06b-e4989d3cbe0b
# ╠═97a4186b-1989-4a20-b60e-b86bb6df9a6a
# ╠═76781e85-099c-4b2c-bb6e-96536523750e
# ╟─74f82d5c-7dbf-4aa6-b021-3a3881d9ad5c
# ╠═1963212f-6dbc-4ca1-aa8c-0f9958a5909d
# ╠═74ccb6be-dd28-45f0-a731-f69a75a17fc3
# ╠═e11c2e4e-f27d-4813-9861-fa2a6f23e523

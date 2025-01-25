### A Pluto.jl notebook ###
# v0.20.3

using Markdown
using InteractiveUtils

# ╔═╡ 30e367ae-a8c8-11ef-2d8b-cd9109046523
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ f15f75d9-ce2b-4048-aa61-45aab8ec5155
begin
	using Plots
	using TaylorModels
	include("plotting.jl")
	include("euler.jl")
	include("tm_integration.jl")
end

# ╔═╡ 8f29564c-a24b-42a3-a27e-748c65b12bee
include("plotting.jl")

# ╔═╡ e1e1828e-4c76-4645-99c7-597eddcd4bbb
"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function vector_field!(du::Vector, u::Vector)
	y, = u

    # The vector field f of the ODEs:
    #   f[1] = y * cos(y)
	du[1] = y * cos(y)
end

# ╔═╡ 834d0782-08ae-4086-bb73-b0de5c3724ef
begin
    # Full TM integration, based on a TMJets bug report:
	# https://github.com/JuliaReach/ReachabilityAnalysis.jl/issues/588
    k = 6   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 10
    NR_REFINEMENTS           = 15
    SCALE                    = 2.0
    TIME_STEP_SIZE = 0.01
    TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
    time_horizon = 1.0     # The finite time horizon
	vars = set_variables("y t", order=k)

	# The variable domains, initial hyperrectangle
    dom::IntervalBox = IntervalBox([
		0.5..2.0,   # x
	])

	# Initial remainder estimate J, a hyperrectangle
	J = IntervalBox(fill(-0.1..0.1, length(dom)))

	PRECONDITIONING::Bool = true
	RECOMPUTE_VECTOR_FIELD::Bool = false

	"k=$k", "NR_REFINEMENTS=$NR_REFINEMENTS", "TIME_STEP_SIZE=$TIME_STEP_SIZE", "time_horizon=$time_horizon", "domains=$dom", "typeof(dom)=$(typeof(dom))" # Pluto cell output
end

# ╔═╡ 262ee75c-7d82-4bd7-a475-85031c80e609
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

# ╔═╡ 4c851641-66f5-4782-a4ad-f1eccdc6b857
begin
	# If true, then plots the boxes computed using the domains that are recomputed every iteration of TM integration.
	# If false, then plots the boxes computed using exclusively the initial domains/hyperrectangle.
	PLOT_DYNAMIC_BOXES::Bool = false

	draw_boxes = PLOT_DYNAMIC_BOXES ? dynamic_boxes : initial_boxes
end

# ╔═╡ 1b7173f1-87b0-4e47-970a-dbe5fbebff25
begin
	# Write output to disk
	filename::String = "TMJets_k$(k)_t$(time_horizon)_dt$(TIME_STEP_SIZE)_R$(NR_REFINEMENTS)"
	open("output/$filename.txt", "w") do io
		for b in draw_boxes
			write(io, "$b\n")
		end
	end
end

# ╔═╡ f17d3b1f-ac1c-413d-b168-74ea2cfcdbb7
begin
	# The number of flowpipes/boxes to keep.
	trunc_max = 50
	trunc_max = length(draw_boxes)

	# Truncate several series, so that you can play around with the number
	# of flowpipes/boxes to display on the plots without having to re-run
	# TM integration.
	time_horizon_truncated = trunc_max * TIME_STEP_SIZE
	step_sizes_truncated = step_sizes[1:trunc_max]
	draw_boxes_truncated = initial_boxes[1:trunc_max]
	known_boxes_truncated::Array{IntervalBox} = zero(draw_boxes_truncated)
end

# ╔═╡ e5b9c19c-14ee-443b-b74e-c988897f1ea5
"""The Forward Euler vector field."""
function vector_field_euler!(du, u, p, t)
	return vector_field!(du, u)
end

# ╔═╡ c520f557-2c51-46c2-984d-3d92e39f6100
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

# ╔═╡ 9bec83c3-8b5d-4a47-93ec-db295178044e
begin
	# Actual plotting
	pltND = plot_boxes_ND(draw_boxes, step_sizes, vars_no_t)

	plot!(pltND, eseries[1])
end

# ╔═╡ Cell order:
# ╠═30e367ae-a8c8-11ef-2d8b-cd9109046523
# ╠═f15f75d9-ce2b-4048-aa61-45aab8ec5155
# ╠═e1e1828e-4c76-4645-99c7-597eddcd4bbb
# ╠═834d0782-08ae-4086-bb73-b0de5c3724ef
# ╠═262ee75c-7d82-4bd7-a475-85031c80e609
# ╠═4c851641-66f5-4782-a4ad-f1eccdc6b857
# ╠═1b7173f1-87b0-4e47-970a-dbe5fbebff25
# ╠═f17d3b1f-ac1c-413d-b168-74ea2cfcdbb7
# ╠═8f29564c-a24b-42a3-a27e-748c65b12bee
# ╠═e5b9c19c-14ee-443b-b74e-c988897f1ea5
# ╠═c520f557-2c51-46c2-984d-3d92e39f6100
# ╠═9bec83c3-8b5d-4a47-93ec-db295178044e

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
	include("ode_shifting.jl")
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
    TIME_STEP_SIZE = 0.02
    TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
    time_horizon = 5.0     # The finite time horizon

	# The variable domains, initial hyperrectangle
    dom::IntervalBox = IntervalBox([
		0.5..2.0,   # x
	])

	# Initial remainder estimate J, a hyperrectangle
	J = IntervalBox(fill(-0.1..0.1, length(dom)))

	PRECONDITIONING::Bool = true
	RECOMPUTE_VECTOR_FIELD::Bool = false
	PLOT_DYNAMIC_BOXES::Bool = false

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

	if PLOT_DYNAMIC_BOXES
		initial_boxes = dynamic_boxes
	end

	initial_sets, step_sizes, initial_boxes, dynamic_boxes # Pluto cell output
end

# ╔═╡ 1b7173f1-87b0-4e47-970a-dbe5fbebff25
begin
	# Write output to disk
	filename::String = "TMJets_k$(k)_t$(time_horizon)_dt$(TIME_STEP_SIZE)_R$(NR_REFINEMENTS)"
	open("output/$filename.txt", "w") do io
		for b in initial_boxes
			write(io, "$b\n")
		end
	end
end

# ╔═╡ f17d3b1f-ac1c-413d-b168-74ea2cfcdbb7
begin
	# Boolean flags to toggle which subplots to show and how.
	USE_LOCAL_HORIZON = false

	# The number of flowpipes/boxes to keep.
	trunc_max = 50
	trunc_max = length(initial_boxes)

	# Truncate several series, so that you can play around with the number
	# of flowpipes/boxes to display on the plots without having to re-run
	# TM integration.
	time_horizon_truncated = trunc_max * TIME_STEP_SIZE
	step_sizes_truncated = step_sizes[1:trunc_max]
	initial_boxes_truncated = initial_boxes[1:trunc_max]
	known_boxes_truncated::Array{IntervalBox} = zero(initial_boxes_truncated)
end

# ╔═╡ 6f7f3780-cca6-4db1-99e7-dd21a67932bc
begin
    vars_no_t = get_variable_names()[1:end-1]
	euler_init_state = [
		1.0  # y
	]
	# Generate the composed plot, which also does euler
	plt_composed, tseries_, vseries_ = plot_vars_against_time(initial_boxes_truncated, known_boxes_truncated,
		time_horizon_truncated, step_sizes_truncated, vars_no_t, USE_LOCAL_HORIZON, odeTMJets!, euler_init_state, euler_step_size = 0.02)

	plt = plot(plot(plt_composed[1]))
	savefig(plt, "output/$filename.png")
	plt # Pluto cell output
end

# ╔═╡ Cell order:
# ╠═30e367ae-a8c8-11ef-2d8b-cd9109046523
# ╠═f15f75d9-ce2b-4048-aa61-45aab8ec5155
# ╠═e1e1828e-4c76-4645-99c7-597eddcd4bbb
# ╠═834d0782-08ae-4086-bb73-b0de5c3724ef
# ╠═262ee75c-7d82-4bd7-a475-85031c80e609
# ╠═1b7173f1-87b0-4e47-970a-dbe5fbebff25
# ╠═f17d3b1f-ac1c-413d-b168-74ea2cfcdbb7
# ╠═8f29564c-a24b-42a3-a27e-748c65b12bee
# ╠═6f7f3780-cca6-4db1-99e7-dd21a67932bc

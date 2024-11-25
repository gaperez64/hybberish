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
	using ProgressBars
	include("plotting.jl")
	include("euler.jl")
	include("tm_integration.jl")
	include("ode_shifting.jl")
end

# ╔═╡ 8f29564c-a24b-42a3-a27e-748c65b12bee
include("plotting.jl")

# ╔═╡ 834d0782-08ae-4086-bb73-b0de5c3724ef
begin
    # Full TM integration, based on a TMJets bug report:
	# https://github.com/JuliaReach/ReachabilityAnalysis.jl/issues/588
    k = 2   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 5
    NR_REFINEMENTS           = 15
    SCALE                    = 2.0
    TIME_STEP_SIZE = 0.01
    TIME_STEP_SIZE_EPS = 2.0e-8  # The minimum time step-size
    time_horizon = 0.02     # The finite time horizon
    vars = set_variables("y t", order=k)
    vars_no_t = get_variable_names()[1:end-1]

    dom::IntervalBox = IntervalBox([
		1.0..2.0,   # x
    	0.0..(TIME_STEP_SIZE)  # t
	])

	# Shift the domains, and replace all variables by shifted variables.
	dom, vars_tms, shift_offsets = shifted_vars_tms(dom)

    # FIXME: Bumping up the order to avoid assertion errors
    # during TM arithmetic.
    old_order = get_order()
    new_order = old_order*2
    vars = set_variables(get_variable_string(), order=new_order)

    # The vector field f of the ODEs:
    #   f[1] = y * cos(y)
    f = [
		vars_tms[1] * cos(vars_tms[1]),
	]
    # Initial remainder estimate J, a hyperrectangle
    J = fill(-0.1..0.1, length(f))

    # TODO: Reset to old order.
    vars = set_variables(get_variable_string(), order=old_order)

	"k=$k", "NR_REFINEMENTS=$NR_REFINEMENTS", "TIME_STEP_SIZE=$TIME_STEP_SIZE", "time_horizon=$time_horizon", dom, vars_tms, f # Pluto cell output
end

# ╔═╡ 1d076713-7065-4e20-bcb3-4e7536f88bce
initial_sets, step_sizes =
	tm_integration(f, dom, k, J, time_horizon,
				   TIME_STEP_SIZE,
				   TIME_STEP_SIZE_EPS,
				   NR_CONTRACTIVENESS_TRIES,
				   NR_REFINEMENTS,
				   SCALE)

# ╔═╡ e0184281-cda1-4133-8876-578d5dfd57b7
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

# ╔═╡ f17d3b1f-ac1c-413d-b168-74ea2cfcdbb7
begin
	include("scripts/clean_intervals_3_3_11.jl")  # Import most recent known boxes

	# Boolean flags to toggle which subplots to show and how.
	SHOW_ONLY_COMPUTED_BOXES = true
	USE_LOCAL_HORIZON = false

	# The number of flowpipes/boxes to keep.
	trunc_max = 182
	trunc_max = length(initial_boxes)

	# Check intersection between computed and known boxes
	# known_boxes = copy(get_clean_boxes())
	known_boxes = map((box) -> zero(box), initial_boxes)

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

# ╔═╡ 1b7173f1-87b0-4e47-970a-dbe5fbebff25
begin
	# Write output to disk
	open("output/TMJets_k$(k)_t$(time_horizon)_dt$(TIME_STEP_SIZE)_R$(NR_REFINEMENTS).txt", "w") do io
		for b in initial_boxes
			write(io, "$b\n")
		end
	end
end

# ╔═╡ b42ee33f-ad37-4745-95b1-5232bb903bbc
map((b) -> diam(b), initial_boxes)

# ╔═╡ c3d1eb22-93e8-482d-8f77-4349fcf9072b
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

# ╔═╡ 6f7f3780-cca6-4db1-99e7-dd21a67932bc
begin
	euler_init_state = [
		0.1  # y
	]
	# Generate the composed plot, which also does euler
	plt_composed, tseries_, vseries_ = plot_vars_against_time(initial_boxes_truncated, known_boxes_truncated, time_horizon_truncated, step_sizes_truncated, vars_no_t, USE_LOCAL_HORIZON, odeTMJets!, euler_init_state, euler_step_size = 0.00001)

	if SHOW_ONLY_COMPUTED_BOXES
		plot(plot(plt_composed[1]), ylims=(-1,3)) # Pluto cell output
	else
		plt_composed # Pluto cell output
	end
end

# ╔═╡ Cell order:
# ╠═30e367ae-a8c8-11ef-2d8b-cd9109046523
# ╠═f15f75d9-ce2b-4048-aa61-45aab8ec5155
# ╠═834d0782-08ae-4086-bb73-b0de5c3724ef
# ╠═1d076713-7065-4e20-bcb3-4e7536f88bce
# ╠═e0184281-cda1-4133-8876-578d5dfd57b7
# ╠═1b7173f1-87b0-4e47-970a-dbe5fbebff25
# ╠═b42ee33f-ad37-4745-95b1-5232bb903bbc
# ╠═f17d3b1f-ac1c-413d-b168-74ea2cfcdbb7
# ╠═c3d1eb22-93e8-482d-8f77-4349fcf9072b
# ╠═8f29564c-a24b-42a3-a27e-748c65b12bee
# ╠═6f7f3780-cca6-4db1-99e7-dd21a67932bc

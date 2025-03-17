### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 4550bf56-a5f4-11ef-2414-599d7eb209c0
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ b9473998-273d-43de-8c3e-4103c63c9822
begin
	using Plots
	include("euler.jl")
end

# ╔═╡ 9dcf3c6b-7e17-4b1a-9cbb-5542664930de
begin
	function euler_plot(ode_func, init_conditions, title)
		# Compute the forward Euler approximation, which yields
		# 	tseries = all time steps where an approximation was computed
		# 	vseries = [
		# 		the x variable approximations,
		# 		the y variable approximations
		# 	]
		euler_step_size = 0.001
		datax_max = 3.1

		plt = plot(xlabel="x", ylabel="y")
		for init in init_conditions
			_, vseries = euler(ode_func, datax_max, euler_step_size, init)
		
			# Display forward Euler results legibly
			datax = map((vec) -> vec[1], vseries)
			datay = map((vec) -> vec[2], vseries)
			plt = plot!(datax, datay, labels="x=$(init[1]), y=$(init[2])", title=title)
		end
		return plt
	end
end

# ╔═╡ 6c5a9593-d507-4eb0-81e8-451b4c9d3625
begin
	# SX: The amount to shift the x-component on the x-axis
	SX = 5.0
	# SY: The amount to shift the y-component on the y-axis
	SY = 2.0
	function ode_lotka_volterra_shifted!(du, u, p, t)
		x, y = u
		du[1] = 1.5 * (x+SX) - (x+SX)*(y+SY)
		du[2] = -3.0 * (y+SY) + (x+SX)*(y+SY)
	end

	# The original lotka volterra ODE, for comparison.
	function ode_lotka_volterra!(du, u, p, t)
		# x' = 1.5x - xy
		# y' = -3y  + xy
	    du[1] = 1.5u[1] - u[1] * u[2]
	    du[2] = -3 * u[2] + u[1] * u[2]
	end
end

# ╔═╡ 23767e5a-59b8-47da-8777-08b59ce02d94
euler_plot(ode_lotka_volterra!, [
	[5.0; 2.0],
	[2.0; 2.0],
	[0.1; 0.1],
	[0.01; 0.01],
], "oridinal ODEs")

# ╔═╡ e4fb0f1f-2b8c-46a8-8121-7be2a8d750f2
begin
	plt = euler_plot(ode_lotka_volterra_shifted!, [
		[10.0; 4.0],
		[5.0; 2.0],
		[0.0; 0.0],
	], "shifted ODEs")
	plot(plt, xlims=(-6, 16), ylims=(-3, 19))
end

# ╔═╡ e14f1ea7-7782-424f-be3a-309402d37253
begin
	plt1 = euler_plot(ode_lotka_volterra!, [
		[10.0; 4.0],
		[5.0; 2.0],
		[0.0; 0.0],
		[15.0; 6.0],
	], "original ODEs")
	plt2 =	euler_plot(ode_lotka_volterra_shifted!, [
		[10.0; 4.0],
		[5.0; 2.0],
		[0.0; 0.0],
		[-5.0; -2.0],
	], "shifted ODEs")
	plot(plt1, plt2, xlims=(-6, 16), ylims=(-3, 19))
end

# ╔═╡ Cell order:
# ╠═4550bf56-a5f4-11ef-2414-599d7eb209c0
# ╠═b9473998-273d-43de-8c3e-4103c63c9822
# ╠═9dcf3c6b-7e17-4b1a-9cbb-5542664930de
# ╠═6c5a9593-d507-4eb0-81e8-451b4c9d3625
# ╠═23767e5a-59b8-47da-8777-08b59ce02d94
# ╠═e4fb0f1f-2b8c-46a8-8121-7be2a8d750f2
# ╠═e14f1ea7-7782-424f-be3a-309402d37253

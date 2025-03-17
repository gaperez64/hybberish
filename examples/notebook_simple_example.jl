### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 11cd7116-0264-11f0-11fd-c3e6362929dc
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ bd5fbea0-a286-48d9-bb0e-75e1f9c1c469
begin

# Re-import 
using TaylorSeries
using Plots
vnames = get_variable_names()

end

# ╔═╡ cfd3f8bc-60dd-47cf-83ba-1f776e3544d6
begin
    include("../src/plotting.jl")
    include("../src/euler.jl")
    include("../src/tm_integration.jl")
end

# ╔═╡ f0f2d653-0905-4a55-abcd-fefd266b8933
md"""
# Setup
"""

# ╔═╡ a66776d2-e62b-4a7e-8483-3044ae84b842
md"""
# Define the Dynamics

We define a constructor function that is used to
1. Taylorize the dynamics as input to the Lie derivative routine.
2. Perform the TM composition in the TM extension of the picard operator.
"""

# ╔═╡ 1c0ff5c8-bec2-4d65-b14c-462ccb2209ee
"""Construct the dynamics.

    This function is used to:
        1. Taylorize the dynamics as input to the Lie derivative routine.
        2. Perform the TM composition in the TM extension of the picard operator.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function f_dot!(du::Vector, u::Vector)
    y, t = u
    du[1] = -y - sin(t) + cos(t)
end


# ╔═╡ c990a5ab-3a3e-49c4-bea1-490d75993164
md"""
# Specify the TM Integration Config
"""

# ╔═╡ c31a35e1-a42c-4c39-a6f8-73f51e8683b5
begin

# The truncation degree / the degree of all polynomials
# that are used during computations.
ord = 10

# Initial state variable bounds and domain.
# y(0) = [1, 1]
# t(0) = [0, 0]
# This is D_0 specifically in the maths.
initial = [
    ("y", interval(1)),
    ("t", interval(0))
]
init_values = IntervalBox([e[2] for e in initial])

# The fixed time step size.
tstep::Float64 = 0.001
# The number of TM integration algo iterations.
nr_iterations::Integer = 120
# Specify time as a finite time horizon.
time_horizon::Float64 = nr_iterations * tstep

# The scale factor for when contractiveness fails.
nr_contractiveness_tries = 1
# The number of refinements to perform at most.
nr_refinements = 5
# Quit refinement early if the improvement a single refinement
# provides falls below this threshold.
refinement_eps = 0.001
# The scale factor with which to widen the initial safe remainder
# estimate when the contractiveness check fails.
scale = 2.0

end

# ╔═╡ f7c13bb2-7c59-40e6-bb0c-86af212f7913
md"""
# Actually run TM Integration
"""

# ╔═╡ 94354bf4-a3ea-4b50-a464-70de8a5cfa39
begin

boxes::Vector{IntervalBox} = []
fboxes::Vector{IntervalBox} = []
boxes, fboxes = tm_integration(
    f_dot!,
    initial,
    ord,
    (-0.1..0.1),
    time_horizon,
    tstep,
    nr_contractiveness_tries,
    nr_refinements,
    SCALE=scale,
    REFINEMENT_EPS=refinement_eps
)

end

# ╔═╡ ef5e4576-5c96-4d47-86d2-fb79f7bfbd47
md"""
# Plotting

**If you want to truncate the series of boxes to inspect a smaller time horizon, see the cell that displays the plots below.**

Note that the solution to the ODEs is "y(t) = cos(t)".
Given that we know the solution y(t) to the ODEs

```
     y(t) = cos(t)
    y'(t) = d(cos(t))/dt = - sin(t)
```

which we verify against the ODEs
```
    y'(t) = -y(t) - sin(t) + cos(t)
          = -(cos(t)) - sin(t) + cos(t)
          = - sin(t)
```
"""

# ╔═╡ a17f50f8-8733-45f9-b86e-13a30afb135d
"""Construct the vector field of the given ODEs.

	@param[out] du The vector field, the right-hand side of the ODEs.
	@param[in]   u The ODE variables to use in construction.
"""
function ode_euler!(du, u, p, t)
	y, = u
    # This function implements the stable version of these ODEs.
    # The unstable versions is: "y - sin(t) - cos(t)"
	du[1] =  -y - sin(t) + cos(t)
end

# ╔═╡ d35ec5a5-891c-4059-994a-c069c899ed33
begin

# Start Forward Euler in the middle of the variable domains.
euler_init_state = Vector(mid(init_values))

# Evaluate Forward Euler.
euler_step_ratio::Integer = 10
euler_step::Float64 = tstep / euler_step_ratio
eseries = euler(ode_euler!, time_horizon, euler_step, euler_init_state)

end

# ╔═╡ cb7f9c5d-e6a1-4785-a05a-b28dc1d6049d
begin

# (Optional) Truncate the vector of boxes, to inspect
# the TM integration results.
nr_boxes = 80
nr_boxes = length(boxes) # TODO: Comment this line if you want to truncate

@assert nr_boxes <= length(boxes)
@assert nr_boxes <= length(fboxes)

truncated_boxes = boxes[1:nr_boxes]
truncated_fboxes = fboxes[1:nr_boxes]
truncated_eseries = [
	series[1:nr_boxes*euler_step_ratio]
	for series in eseries
]


# Generate the plots
pltND1 = plot_boxes_ND(truncated_boxes, vnames, sgtitle="vals", legend=true)
pltND2 = plot_boxes_ND(truncated_fboxes, vnames, sgtitle="fpipe", legend=true)
plot!(pltND1, truncated_eseries[1], label="Stable ODE Forward Euler")
plot!(pltND1, cos, label="cos(t)") # The ODE solution is "y(t) = cos(t)"
plot!(pltND2, truncated_eseries[1], label="Stable ODE Forward Euler")
plot!(pltND2, cos, label="cos(t)") # The ODE solution is "y(t) = cos(t)"


# Make Pluto notebook display the plots
pltND = plot(pltND1, pltND2)

end

# ╔═╡ f40f83aa-c9ee-4308-a2b9-7f7cbc91bb43
md"""
# Calling the notebook in the CLI

In the CLI, call `julia -i xxx.jl` where you replace "xxx" by the notebook name. The `-i` flag ensures the generated plot is actually displayed.

For example:

`julia -i notebook_simple_example.jl`

"""

# ╔═╡ 4e4db550-a643-465e-90b2-7ed54a440446
display(pltND)

# ╔═╡ Cell order:
# ╟─f0f2d653-0905-4a55-abcd-fefd266b8933
# ╠═11cd7116-0264-11f0-11fd-c3e6362929dc
# ╠═cfd3f8bc-60dd-47cf-83ba-1f776e3544d6
# ╟─a66776d2-e62b-4a7e-8483-3044ae84b842
# ╠═1c0ff5c8-bec2-4d65-b14c-462ccb2209ee
# ╟─c990a5ab-3a3e-49c4-bea1-490d75993164
# ╠═c31a35e1-a42c-4c39-a6f8-73f51e8683b5
# ╟─f7c13bb2-7c59-40e6-bb0c-86af212f7913
# ╠═94354bf4-a3ea-4b50-a464-70de8a5cfa39
# ╟─ef5e4576-5c96-4d47-86d2-fb79f7bfbd47
# ╠═a17f50f8-8733-45f9-b86e-13a30afb135d
# ╠═d35ec5a5-891c-4059-994a-c069c899ed33
# ╠═bd5fbea0-a286-48d9-bb0e-75e1f9c1c469
# ╠═cb7f9c5d-e6a1-4785-a05a-b28dc1d6049d
# ╟─f40f83aa-c9ee-4308-a2b9-7f7cbc91bb43
# ╠═4e4db550-a643-465e-90b2-7ed54a440446

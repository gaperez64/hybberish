### A Pluto.jl notebook ###
# v0.20.1

using Markdown
using InteractiveUtils

# ╔═╡ 59be0dc5-fe79-465d-918e-da5778cddc54
begin
    import Pkg
    # I used Base.current_project() because you mentioned that your notebook is located in the same folder of the Project.toml
    Pkg.activate()
end

# ╔═╡ 75e8220e-f7bb-4e1a-8c51-0e0860593c0f
begin
	using TaylorSeries  # IntervalArithmetic

	println(pathof(TaylorSeries))
	
	# Assumes the last variable is t
	function tay_poly(f, k)
	    vars = get_variables()
	    t = vars[end]
	    # Let's pad f with a 1 at the end for t
	    fp1 = copy(f)
	    push!(fp1, 1)
	    # Also, prepare a valuation vector with t=0
	    # FIXME: We're cheating to obtain a zero with same order as the other
	    # variables
	    val0 = copy(vars)
	    val0[end] = t - t
	    # Prepare a first lie derivative and the result
	    g = copy(vars)
	    deleteat!(g, length(vars))
	    # Start a vector function for the result
	    res = copy(g)
	    
	    for i = 1:k
	        println("g = $g")
	        g = TaylorSeries.jacobian(g, vars) * fp1
	        term = map((h) -> evaluate(h, val0) * t^i * (1 / factorial(i)), g)
	        res += term
	    end
	    
	    return res
	end
end

# ╔═╡ 0c67b6fe-9163-11ef-019e-333f0ff3dc9d
md"""
## Introduction
Text can be added in Markdown format. Rendering is done with the "play" button, CMD+S, CMD+Enter, ...
"""

# ╔═╡ 7066cb63-a64d-40e8-8a0e-684278b24e98
# ╠═╡ show_logs = false
# add some code
println("Hello world")

# ╔═╡ ef7165a0-b7ca-4836-b8d3-7d557555efb5
md"""
Be sure to check out the "Live Docs" and "Status" buttons as well.
"""

# ╔═╡ e6f4bcd2-7b76-4b38-b52b-bbf5ce52ea51
md"""## The actual work begins ...
Below a function called `tay_poly` is defined. It takes a function $f$, and an integer $k$.

Make sure you use the version of `TaylorSeries` that supports non-square jacobians.
"""

# ╔═╡ 3de35676-e988-4ad4-8b91-fac9507aebdf


# ╔═╡ b8eafb68-c0f4-4824-a935-71e50f32dbb8
md"""
We now try a small example.
"""

# ╔═╡ dcefa9d1-0d86-4e78-a848-b702aff63fcc
# ╠═╡ show_logs = false
begin
	vars = set_variables("x y t", order=4)
	f = [1 + vars[2], -vars[1]^2]
	p = tay_poly(f, 4)
end

# ╔═╡ e5c18132-3337-4a91-8772-8a1e7d64f486
println("poly from new TM = $p")

# ╔═╡ df802817-ccb1-4e8e-bc65-03d8b4a488bc
md"""
Done.
"""

# ╔═╡ Cell order:
# ╠═0c67b6fe-9163-11ef-019e-333f0ff3dc9d
# ╠═7066cb63-a64d-40e8-8a0e-684278b24e98
# ╠═ef7165a0-b7ca-4836-b8d3-7d557555efb5
# ╟─e6f4bcd2-7b76-4b38-b52b-bbf5ce52ea51
# ╠═59be0dc5-fe79-465d-918e-da5778cddc54
# ╠═3de35676-e988-4ad4-8b91-fac9507aebdf
# ╠═75e8220e-f7bb-4e1a-8c51-0e0860593c0f
# ╠═b8eafb68-c0f4-4824-a935-71e50f32dbb8
# ╠═dcefa9d1-0d86-4e78-a848-b702aff63fcc
# ╠═e5c18132-3337-4a91-8772-8a1e7d64f486
# ╠═df802817-ccb1-4e8e-bc65-03d8b4a488bc

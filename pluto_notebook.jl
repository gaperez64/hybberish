### A Pluto.jl notebook ###
# v0.20.1

using Markdown
using InteractiveUtils

# ╔═╡ 75e8220e-f7bb-4e1a-8c51-0e0860593c0f
begin
	using TaylorSeries  # IntervalArithmetic
	
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

# ╔═╡ bc458c4d-17ed-45a3-8301-2040e414cf76


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
	println("poly from new TM = $p")
end

# ╔═╡ df802817-ccb1-4e8e-bc65-03d8b4a488bc
md"""
Done.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
TaylorSeries = "6aa5eb33-94cf-58f4-a9d0-e4b2c4fc25ea"

[compat]
TaylorSeries = "~0.18.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.1"
manifest_format = "2.0"
project_hash = "e0f210c458782bc103fc2c8c91a33988e41db555"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TaylorSeries]]
deps = ["LinearAlgebra", "Markdown", "Requires", "SparseArrays"]
git-tree-sha1 = "bb212ead98022455eed514cff0adfa5de8645258"
uuid = "6aa5eb33-94cf-58f4-a9d0-e4b2c4fc25ea"
version = "0.18.1"

    [deps.TaylorSeries.extensions]
    TaylorSeriesIAExt = "IntervalArithmetic"
    TaylorSeriesJLD2Ext = "JLD2"
    TaylorSeriesRATExt = "RecursiveArrayTools"
    TaylorSeriesSAExt = "StaticArrays"

    [deps.TaylorSeries.weakdeps]
    IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
    JLD2 = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
    RecursiveArrayTools = "731186ca-8d62-57ce-b412-fbd966d074cd"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"
"""

# ╔═╡ Cell order:
# ╠═0c67b6fe-9163-11ef-019e-333f0ff3dc9d
# ╠═bc458c4d-17ed-45a3-8301-2040e414cf76
# ╠═7066cb63-a64d-40e8-8a0e-684278b24e98
# ╠═ef7165a0-b7ca-4836-b8d3-7d557555efb5
# ╟─e6f4bcd2-7b76-4b38-b52b-bbf5ce52ea51
# ╠═75e8220e-f7bb-4e1a-8c51-0e0860593c0f
# ╠═b8eafb68-c0f4-4824-a935-71e50f32dbb8
# ╠═dcefa9d1-0d86-4e78-a848-b702aff63fcc
# ╠═df802817-ccb1-4e8e-bc65-03d8b4a488bc
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

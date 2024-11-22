using YAML.jl

"""Assumes that all activation functions are ReLU"""
function load_dnn_fun(yml_filename)
    data = YAML.load_file(yml_filename)
    offsets = [data["offsets"][key] for key in 1:length(data["offsets"])]
    weights = [data["weights"][key] for key in 1:length(data["weights"])]
    zs = [zeros(length(offsets[i])) for i in eachindex(offsets)]
    
    funs = [x -> max.(weights[i] * x + offsets[i], zs[i])
            for i in eachindex(offsets)]

    return reduce(|>, funs)
end


if abspath(PROGRAM_FILE) == @__FILE__
    # Example 3.3.6
    k = 3   # The TM arithmetic and truncation order
    NR_CONTRACTIVENESS_TRIES = 5
    NR_REFINEMENTS           = 1
    SCALE                    = 2.0
    vars = set_variables("x y t", order=k)
    # The vector field f of the ODEs:
    #   f[1] = 1 + y
    #   f[2] = -x^2
    f = [1 + vars[2],  # x
         -vars[1]^2]   # y
    # The polynomial approx:
    #   p[1] = x + t + yt
    #   p[2] = y - (x^2)t - xt^2 - (1/3)t^3
    p = tay_poly(f, k)  # FIXME: We should be using a general fun!
    domain = IntervalBox([-1..1,      # x
    		      -0.5..0.5,  # y
    		      0..0.02])   # t
    # Initial remainder estimate J, a hyperrectangle
    J = fill(-0.1..0.1, length(f))
    # Let's get that safe remainder now!
    I = tay_model_error(f, p, domain, k, J,
                        NR_CONTRACTIVENESS_TRIES,
                        NR_REFINEMENTS,
                        SCALE)
    println("safe remainders = $I")
end

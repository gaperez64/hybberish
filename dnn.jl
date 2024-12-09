using YAML

"""Assumes that all activation functions are ReLU"""
function load_dnn_fun(yml_filename)
    data = YAML.load_file(yml_filename)
    offsets = [data["offsets"][key] for key in 1:length(data["offsets"])]
    weights = [stack(data["weights"][key], dims=1)
               for key in 1:length(data["weights"])]
    zs = [zeros(length(offsets[key])) for key in 1:length(data["offsets"])]
    funs = [i < length(offsets) ?
            x -> max.(weights[i] * x + offsets[i], zs[i]) :
            x -> (weights[i] * x + offsets[i])
            for i in eachindex(offsets)]
    return reduce(∘, reverse(funs))
end


if abspath(PROGRAM_FILE) == @__FILE__
    dnn_filename = ARGS[1]
    println("Loading ReLU DNN stored in $dnn_filename")
    f = load_dnn_fun(dnn_filename)
    vals = map(x -> parse(Float64, x), ARGS[2:end])
    println("Arguments for the DNN function: $vals")
    res = f(vals)
    println("Result: $res")
end

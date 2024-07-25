struct SpeedyOutput{FN, B} # I will make this an `AbstractOutputWriter` at some point but this will do for now
              filename :: FN
    overwrite_existing :: B

    function SpeedyOutput(filename::FN; overwrite_existing::B = true) where {FN, B}
        jldopen(filename, ifelse(overwrite_existing, "w+", "w")) do file
            file["filename"] = filename
        end

        return new{FN, B}(filename, overwrite_existing)
    end
end

function (save::SpeedyOutput)(simulation)
    model = simulation.model

    t = time(simulation)
    iter = model.clock.iteration

    jldopen(save.filename, "a+") do file

        file["timeseries/t/$iter"] = t

        for (i, name) in enumerate(keys(model.fields))
            file["timeseries/$name/$iter"] = model.field_values[i]
        end
    end

    return nothing
end

function load_output(save::SpeedyOutput)
    names = tuple()
    
    jldopen(save.filename) do f 
        names = Symbol.(keys(f["timeseries"]))
    end

    sort_order = sortperm(load_output(save, "t"))

    return NamedTuple{tuple(names...)}([load_output(save, name) for name in names])
end
    
function load_output(save::SpeedyOutput, name)
    sort_order = output_order(save)

    return Float64.(values(load(save.filename, "timeseries/"*string(name))))[sort_order]
end

load_output(save::SpeedyOutput, names::Tuple) = 
    NamedTuple{name}([load_output(save, name) for name in names])

output_order(save) = sortperm(Float64.(values(load(save.filename, "timeseries/t"))))
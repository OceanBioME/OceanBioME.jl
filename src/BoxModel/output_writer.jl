struct SpeedyOutput{FN, B} # I will make this an `AbstractOutputWriter` at some point but this will do for now
              filename :: FN
    overwrite_existing :: B

    function SpeedyOutput(filename::FN; overwrite_existing::B = true) where {FN, B}
        return new{FN, B}(filename, overwrite_existing)
    end
end

function (save::SpeedyOutput)(simulation)
    model = simulation.model

    t = time(simulation)
    iter = model.clock.iteration

    file = jldopen(save.filename, ifelse(save.overwrite_existing, "w+", "w"))

    file["timeseries/t/$iter"] = t

    for (i, name) in enumerate(keys(model.fields))
        file["timeseries/$name/$iter"] = model.field_values[i]
    end

    close(file)

    return nothing
end
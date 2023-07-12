using OceanBioME, Test, DataDeps, JLD2

const minute = minutes = 60
const hour = hours = 60 * minutes
const day = days = hours * 24 

function load_boxmodel(path, vars)
    file = jldopen(path)
    times = parse.(Float64, keys(file["values"]))

    timeseries = NamedTuple{vars}(ntuple(t -> zeros(length(times)), length(vars)))

    for (idx, time) in enumerate(times)
        values = file["values/$time"]
        for tracer in vars
            getproperty(timeseries, tracer)[idx] = values[tracer]
        end
    end

    close(file)

    return times, timeseries
end

function visulise_regression(path, vars)
    times, timeseries = load_boxmodel(path, vars)

    fig = Figure(resolution = (1600, 1600), fontsize = 24)

    axs = []

    len = floor(Int, length(vars) / 2)

    tick = 1
    for (idx, tracer) in enumerate(vars)
        push!(axs, Axis(fig[floor(Int, (idx - 1) / len), tick], ylabel = "$tracer", xlabel = "Days"))
        lines!(axs[end], times / days, timeseries[tracer], linewidth = 3)

        tick = ifelse(tick == 1, 2, 1)
    end

    return fig
end

include("LOBSTER.jl")
include("NPZD.jl")
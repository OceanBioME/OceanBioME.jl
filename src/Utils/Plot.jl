module Plot
using Plots, Oceananigans, JLD2, Statistics
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

load_tracers(sim::Simulation, save_name=:profiles) = load_tracers(sim.output_writers[save_name].filepath)

mutable struct model_results{T, X, Y, Z, t, R}
    tracers :: T
    x :: X
    y :: Y
    z :: Z
    t :: t
    results :: R
end

mutable struct particle_results{P, t, R}
    properties :: P
    t :: t
    results :: R
end

function load_tracers(path)
    file_profiles = jldopen(path)

    tracers=keys(file_profiles["timeseries"])
    deleteat!(tracers, findall(x->x in ["u", "v", "w", "t", "b"], tracers))
    x=file_profiles["grid/xᶜᵃᵃ"][4:end-3]
    y=file_profiles["grid/yᵃᶜᵃ"][4:end-3]
    z=file_profiles["grid/zᵃᵃᶜ"][4:end-3]

    iterations = parse.(Int, keys(file_profiles["timeseries/t"]))

    times = [file_profiles["timeseries/t/$iter"] for iter in iterations]

    results = zeros(length(tracers), length(x), length(y), length(z), length(iterations))
    times = zeros(length(iterations))
    for (i, iter) in enumerate(iterations)
        times[i] = file_profiles["timeseries/t/$iter"]
        for (j, tracer) in enumerate(tracers)
            results[j, :, :, :, i] .= file_profiles["timeseries/$tracer/$iter"]#permutedims(repeat(file_profiles["timeseries/$tracer/$iter"],1,1,1,1,1), (5, 1, 2, 3, 4))
        end
    end
    return model_results(tracers, x, y, z, times, results)
end

load_particles(sim::Simulation, save_name=:particles) = load_particles(sim.output_writers[save_name].filepath)

function load_particles(path)
    file_profiles = jldopen(path)

    iterations = parse.(Int, keys(file_profiles["timeseries/t"]))
    tracers=keys(file_profiles["timeseries/particles/$(iterations[1])"])

    times = [file_profiles["timeseries/t/$iter"] for iter in iterations]

    results = zeros(length(tracers), length(getproperty(file_profiles["timeseries/particles/$(iterations[1])"], tracers[1])), length(iterations))
    times = zeros(length(iterations))
    for (i, iter) in enumerate(iterations)
        times[i] = file_profiles["timeseries/t/$iter"]
        for (j, tracer) in enumerate(tracers)
            results[j, :, i] .= getproperty(file_profiles["timeseries/particles/$iter"], tracer)
        end
    end
    return particle_results(tracers, times, results)
end

function profiles(results::model_results, fs=4, xlabel="time (days)", ylabel="z (m)")
    plts=[]
    for (j, tracer) in enumerate(results.tracers)
        push!(plts, heatmap(results.t/(1day),results.z,mean(results.results[j, :, :, :, :], dims=(1, 2))[1, 1, :, :], titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel, ylabel=ylabel, title=tracer))
    end
    plot(plts...)
end

function particles(results::particle_results, fs=4, xlabel="time (days)")
    plts=[]
    for (j, tracer) in enumerate(results.properties)
        push!(plts, plot(results.t/(1day), results.results[j, :, :]', titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel, ylabel=tracer))
    end
    plot(plts...)
end
end
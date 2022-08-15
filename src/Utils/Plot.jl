module Plot
using Plots, Oceananigans, JLD2, Statistics, NetCDF
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

struct NetCDFResults
    path
    tracers
    tstart
end

function load_tracers(path::String, tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD), tstart=0)
    if path[end-4:end] == ".jld2"
        file_profiles = jldopen(path)
    elseif path[end-2:end] == ".nc"
        file_profiles = NetCDFResults(path, tracers, tstart)
    else
        throw(ArgumentError("$path is not a valid path"))
    end

    return load_tracers(file_profiles)
end

function load_tracers(file_profiles::JLD2.JLDFile)
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
            results[j, :, :, :, i] .= file_profiles["timeseries/$tracer/$iter"]
        end
    end
    return model_results(tracers, x, y, z, times, results)
end

function load_tracers(file::NetCDFResults)
    x = ncread("nonlin_grid.nc", "xC")
    y = ncread("nonlin_grid.nc", "yC")
    z = ncread("nonlin_grid.nc", "zC")
    t = ncread("nonlin_grid.nc", "time")[file.tstart:end]

    results = zeros(length.([file.tracers, x, y, z, t])...)

    for (i, tracer) in enumerate(file.tracers)
        if !(tracer in ["u", "v", "w", "t", "b"])
            results[i, :, :, :, :] =ncread(file.path, "$tracer")[:, :, :, file.tstart:end]
        end
    end

    return model_results(file.tracers, x, y, z, t, results)
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

function profiles(results::model_results; fs=4, xlabel="time (days)", ylabel="z (m)", sediment=("Nᵣ", "Nᵣᵣ", "Nᵣₑ"))
    plts=[]
    for (j, tracer) in enumerate(results.tracers)
        if !(tracer in sediment)
            push!(plts, heatmap(results.t/(1day),results.z,mean(results.results[j, :, :, :, :], dims=(1, 2))[1, 1, :, :], titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel, ylabel=ylabel, title=tracer))
        else
            push!(plts, plot(results.t/(1day),mean(results.results[j, :, :, 1, :], dims=(1, 2))[1, 1, :], titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel, ylabel=ylabel, title=tracer, legend=false))
        end
    end
    plot(plts...)
    return plts
end

function particles(results::particle_results, fs=4, xlabel="time (days)")
    plts=[]
    for (j, tracer) in enumerate(results.properties)
        push!(plts, plot(results.t/(1day), results.results[j, :, :]', titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel, ylabel=tracer, legend=false))
    end
    plot(plts...)
    return plts
end
end
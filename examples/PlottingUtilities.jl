using Plots, Oceananigans, JLD2, Statistics, NetCDF, GLMakie
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using OceanBioME
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

function load_tracers(path::String, tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD), tstart=1)
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
    x = ncread(file.path, "xC")
    y = ncread(file.path, "yC")
    z = ncread(file.path, "zC")
    t = ncread(file.path, "time")[file.tstart:end]

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

function profiles(results::model_results; fs=4, xlabel="time (days)", ylabel="z (m)", sediment=(:Nᵣ, :Nᵣᵣ, :Nᵣₑ))
    plts=[]
    for (j, tracer) in enumerate(results.tracers)
        if (!(tracer in sediment) && !(tracer in ["$t" for t in sediment]))
            push!(plts, Plots.heatmap(results.t/(1day),results.z,mean(results.results[j, :, :, :, :], dims=(1, 2))[1, 1, :, :], titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel, ylabel=ylabel, title=tracer))
        else
            push!(plts, Plots.plot(results.t/(1day),mean(results.results[j, :, :, 1, :], dims=(1, 2))[1, 1, :], titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel, ylabel=ylabel, title=tracer, legend=false))
        end
    end
    Plots.plot(plts...)
    return plts
end

function plot_PISCES_section(results::model_results; fs=4, xlabel="time (days)", ylabel="z (m)", prefix="PISCES_example")
    plts1=[]#Phytoplankton
    plts2=[]#zooplankton and particles
    plts3=[]#Nutrients
    plts4=[]#Chemistry
    for i=1:length(results.tracers)
        push!(i<=7 ? plts1 : (i<=14 ? plts2 : (i<=19 ? plts3 : plts4)), Plots.heatmap(results.t./days, results.z, results.results[i, 1, 1, :, :], title=results.tracers[i]))
    end
    Plots.plot(plts1...)
    savefig(prefix*"_phyto.pdf")
    Plots.plot(plts2...)
    savefig(prefix*"_zoo_particles.pdf")
    Plots.plot(plts3...)
    savefig(prefix*"_nutrients.pdf")
    Plots.plot(plts4...)
    savefig(prefix*"_chem.pdf")
    return plts1, plts2, plts3, plts4
end

function particles(results::particle_results, fs=4, xlabel="time (days)")
    plts=[]
    for (j, tracer) in enumerate(results.properties)
        push!(plts, Plots.plot(results.t/(1day), results.results[j, :, :]', titlefontsize=fs, guidefontsize=fs, tickfontsize=fs, legendfontsize=fs, xlabel=xlabel, ylabel=tracer, legend=false))
    end
    Plots.plot(plts...)
    return plts
end

function plot_particles(path)
    @warn "GLMakie has a memory leak error if you call save lots of times (and various other things)"

    res_kelp = load_particles(path)

    xs, ys = res_kelp.t/(1day), res_kelp.results[3,:,1]

    f = Figure(backgroundcolor = RGBf(1, 1, 1), resolution = (1920, 1050))
    gInput = f[1, 1] = GridLayout()
    gProperties = f[1, 2] = GridLayout()
    gOutput = f[1, 3] = GridLayout()

    ax1, hm1 = GLMakie.heatmap(gInput[1, 1], xs, ys, res_kelp.results[16,:,:]')
    ax1.title = "NO₃"
    ax1.xticklabelsvisible= false
    cb1 = Colorbar(gInput[1, 1:2], hm1, label = "mmol N/m³")

    ax2, hm2 = GLMakie.heatmap(gInput[2, 1], xs, ys, res_kelp.results[17,:,:]')
    ax2.title = "NH₄"
    ax2.ylabel = "depth (m)"
    ax2.xticklabelsvisible= false
    cb2 = Colorbar(gInput[2, 1:2], hm2, label = "mmol N/m³")

    ax3, hm3 = GLMakie.heatmap(gInput[3, 1], xs, ys, res_kelp.results[18,:,:]')
    ax3.title = "PAR"
    ax3.xlabel = "time (day)"
    cb3 = Colorbar(gInput[3, 1:2], hm3, label = "einstein/m²/day")

    ax1, hm1 = GLMakie.heatmap(gProperties[1, 1], xs, ys, res_kelp.results[7,:,:]')
    ax1.title = "Area"
    ax1.xticklabelsvisible= false
    cb1 = Colorbar(gProperties[1, 1:2], hm1, label = "dm²")

    ax2, hm2 = GLMakie.heatmap(gProperties[2, 1], xs, ys, ((res_kelp.results[8,:,:].+SLatissima.defaults.N_struct).*SLatissima.defaults.K_A .*res_kelp.results[7, :, :])')
    ax2.title = "Total Nitrogen (structural + reserve)"
    ax2.ylabel = "depth (m)"
    ax2.xticklabelsvisible= false
    cb2 = Colorbar(gProperties[2, 1:2], hm2, label = "gN/frond")

    ax3, hm3 = GLMakie.heatmap(gProperties[3, 1], xs, ys, ((res_kelp.results[9,:,:].+SLatissima.defaults.C_struct).*SLatissima.defaults.K_A .*res_kelp.results[7, :, :])')
    ax3.title = "Total Carbon (structural + reserve)"
    ax3.xlabel = "time (day)"
    cb3 = Colorbar(gProperties[3, 1:2], hm3, label = "gC/frond")

    ax1, hm1 = GLMakie.heatmap(gOutput[1, 1], xs, ys, res_kelp.results[10,:,:]')
    ax1.title = "NO₃ uptake"
    cb1 = Colorbar(gOutput[1, 1:2], hm1, label = "mmol N/s")
    ax1.xticklabelsvisible= false

    ax2, hm2 = GLMakie.heatmap(gOutput[2, 1], xs, ys, res_kelp.results[11,:,:]')
    ax2.title = "NH₄ uptake"
    ax2.xticklabelsvisible= false
    cb2 = Colorbar(gOutput[2, 1:2], hm2, label = "mmol N/s")

    ax3, hm3 = GLMakie.heatmap(gOutput[3, 1], xs, ys, res_kelp.results[12,:,:]')
    ax3.title = "Primary production (photosynthesis - respiration)"
    ax3.xticklabelsvisible= false
    ax3.ylabel = "depth (m)"
    cb3 = Colorbar(gOutput[3, 1:2], hm3, label = "mmol C/s")

    ax4, hm4 = GLMakie.heatmap(gOutput[4, 1], xs, ys, res_kelp.results[14,:,:]')
    ax4.title = "Exudation (DOC output)"
    ax4.xticklabelsvisible= false
    cb4 = Colorbar(gOutput[4, 1:2], hm4, label = "mmol C/s")

    ax5, hm5 = GLMakie.heatmap(gOutput[5, 1], xs, ys, res_kelp.results[15,:,:]')
    ax5.title = "Frond errosion (POM output)"
    ax5.xlabel = "time (day)"
    cb5 = Colorbar(gOutput[5, 1:2], hm5, label = "mmol N/s")

    save("$path.png", f)
    return f
end
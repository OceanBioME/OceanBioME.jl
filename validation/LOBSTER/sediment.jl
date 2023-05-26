using OceanBioME, Oceananigans, Test, JLD2
using OceanBioME.Sediments: SimpleMultiG
using Oceananigans.Units

year = years = 365days # just for these idealised cases

using CairoMakie

function test_flat_sediment(architecture; timestepper = :QuasiAdamsBashforth2)
    grid = RectilinearGrid(architecture; size=(3, 3, 3), extent=(10, 10, 200))

    sediment_model = SimpleMultiG(grid)

    biogeochemistry = LOBSTER(; grid, 
                                carbonates = true, oxygen = true, variable_redfield = true, 
                                open_bottom = true, 
                                surface_phytosynthetically_active_radiation = (x, y, t) -> 80,
                                sediment_model)

    model = NonhydrostaticModel(;grid, biogeochemistry, 
                                 boundary_conditions = (DIC = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂)), 
                                                        O₂ = FieldBoundaryConditions(top = GasExchange(; gas = :O₂))),
                                 tracers = (:T, :S),
                                 closure = ScalarDiffusivity(ν = 10 ^ -3, κ = 10 ^ -3),
                                 timestepper)

    set!(model.biogeochemistry.sediment_model.fields.N_fast, 0.0230)
    set!(model.biogeochemistry.sediment_model.fields.N_slow, 0.0807)

    set!(model.biogeochemistry.sediment_model.fields.C_fast, 0.5893)
    set!(model.biogeochemistry.sediment_model.fields.C_slow, 0.1677)

    set!(model, P = 0.4686, Z = 0.5363, 
            NO₃ = 2.3103, NH₄ = 0.0010, 
            DIC = 2106.9, Alk = 2408.9, 
            O₂ = 258.92, 
            DOC = 5.3390, DON = 0.8115,
            sPON = 0.2299, sPOC = 1.5080,
            bPON = 0.0103, bPOC = 0.0781)

    simulation = Simulation(model, Δt = 800.0, stop_time = 10years)

    simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                           filename = "sediment_test_tracers_rk3.jld2",
                                                           schedule = TimeInterval(10minutes),
                                                           overwrite_existing = true)

    simulation.output_writers[:sediment] = JLD2OutputWriter(model, model.biogeochemistry.sediment_model.fields,
                                                            indices = (:, :, 1),
                                                            filename = "sediment_test_sediment_rk3.jld2",
                                                            schedule = TimeInterval(10minutes),
                                                            overwrite_existing = true)

    @inline progress(simulation) = @info "Time: $(prettytime(simulation.model.clock.time)), Iteration: $(simulation.model.clock.iteration), Walltime: $(prettytime(simulation.run_wall_time))"
    
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON))
    simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())
    
    plankton_redfield = model.biogeochemistry.phytoplankton_redfield
    scale_negative_carbon_tracers = ScaleNegativeTracers(; model, tracers = (:P, :Z, :DOC, :sPOC, :bPOC, :DIC), 
                                                         scalefactors = (P = plankton_redfield, 
                                                                         Z = plankton_redfield, 
                                                                         DOC = 1, sPOC = 1, bPOC = 1, DIC = 1))
    simulation.callbacks[:neg_carbon] = Callback(scale_negative_carbon_tracers; callsite = UpdateStateCallsite())

    return simulation, model
end

function load_timeseries(; simulation, tracer_path = simulation.output_writers[:tracers].filepath, sediment_path = simulation.output_writers[:sediment].filepath, end_idx = nothing)
    tracer_timeseries = NamedTuple()

    file = jldopen(tracer_path)

    iterations = keys(file["timeseries/t"])

    if isnothing(end_idx) end_idx = length(iterations) end

    for tracer in Symbol.(keys(file["timeseries"]))
        if !(tracer in [:T, :S, :t])
            @info tracer
            data = zeros(end_idx, 3)
            valid_its = keys(file["timeseries/$tracer"])
            for (i, it) in enumerate(iterations[1:end_idx])
                if it in valid_its
                    data[i, :] = file["timeseries/$tracer/$it"][2, 2, 1:3]
                else
                    data[i, :] .= NaN
                end
            end
            tracer_timeseries = merge(tracer_timeseries, NamedTuple{(tracer, )}((data, )))
        end
    end

    zn = znodes(Center, file["serialized/grid"])[1:file["serialized/grid"].Nz]

    close(file)

    sediment_timeseries = NamedTuple()

    file = jldopen(sediment_path)

    iterations = keys(file["timeseries/t"])

    for tracer in Symbol.(keys(file["timeseries"]))
        if !(tracer in [:T, :S, :t])
            @info tracer
            data = zeros(end_idx)
            for (i, it) in enumerate(iterations[1:end_idx])
                data[i] = file["timeseries/$tracer/$it"][2, 2, 1]
            end
            sediment_timeseries = merge(sediment_timeseries, NamedTuple{(tracer, )}((data, )))
        end
    end

    times = zeros(end_idx)
    for (i, it) in enumerate(iterations[1:end_idx])
        times[i] = file["timeseries/t/$it"]
    end

    close(file)

    return tracer_timeseries, sediment_timeseries, times, zn
end

function plot_timeseries(simulation, tracer_timeseries, sediment_timeseries; times = FieldTimeSeries(simulation.output_writers[:tracers].filepath, "P").times, zn = znodes(Center, simulation.model.grid)[1:simulation.model.grid.Nz], start_idx = 1)
    fig = Figure(resolution = (1600, 1600))

    # since they seems to be wrong ...
    times = [1:size(tracer_timeseries.P)[1];] * 10minutes/day

    for (idx, (name, values)) in enumerate(pairs(tracer_timeseries))
        @info "$name"
        ax = Axis(fig[idx % 5, ceil(Int, idx / 5) * 2 - 1], title = "$name", xlabel = "Time (days)", ylabel = "Depth (m)")
        hm = heatmap!(ax, times[start_idx:end], zn, values[start_idx:end, :], colorrange = (minimum(values[start_idx:end, :][isfinite.(values[start_idx:end, :])]), maximum(values[start_idx:end, :][isfinite.(values[start_idx:end, :])])))
        Colorbar(fig[idx % 5, ceil(Int, idx / 5) * 2], hm, label = "$name Concentration (mmol / m³)")
    end

    for (idx, (name, values)) in enumerate(pairs(sediment_timeseries))
        @info "$name"
        ax = Axis(fig[(idx + length(tracer_timeseries)) % 5, ceil(Int, (idx + length(tracer_timeseries)) / 5) * 2 - 1], title = "$name", xlabel = "Time (days)", ylabel = "Concentration (mmol / m²)")
        hm = lines!(ax, times[start_idx:end], values[start_idx:end])
    end

    return fig
end
#=
fig = Figure(resolution = (1600, 1000))
volume(fig[1, 1], xnodes(Face, u.grid)[1:u.grid.Nx], ynodes(Center, u.grid)[1:u.grid.Ny], znodes(Center, u.grid)[1:u.grid.Nz], u_plt, colorrange = (-uₘ, uₘ))
GLMakie.record(fig, "ovs_fig.mp4", 1:length(u.times)) do i
    n[] = i
    msg = string("Plotting frame ", i, " of ", length(u.times))
    print(msg * " \r")
    #ax.title = "t=$(prettytime(u.times[i]))"
end=#

       fig = Figure(resolution = (1600, 1600))
       l = Label(fig[1, 1:4], "t = 0.0")
       v1 = volume(fig[2, 1], xnodes(grid, Face(), Center(), Center())[1:grid.Nx], ynodes(grid, Center(), Center(), Center())[1:grid.Ny], znodes(grid, Center(), Center(), Center())[1:grid.Nz], u_plt, colorrange = (-uₘ, uₘ), colormap=:vik); Colorbar(fig[2, 2], v1.plot, label = "u (m/s)")
       v2 = volume(fig[3, 1], xnodes(grid, Center(), Center(), Center())[1:grid.Nx], ynodes(grid, Center(), Center(), Center())[1:grid.Ny], znodes(grid, Center(), Center(), Center())[1:grid.Nz], NO₃_plt, colorrange = (Nₘᵢ, Nₘₐ), colormap = Reverse(:bamako), algorithm = :mip)
       Colorbar(fig[3, 2], v2.plot, label = "NO₃ (mmol N / m³)")
       v3 = volume(fig[2, 3], xnodes(grid, Center(), Center(), Center())[1:grid.Nx], ynodes(grid, Center(), Center(), Center())[1:grid.Ny], znodes(grid, Center(), Center(), Center())[1:grid.Nz], P_plt, colorrange = (Pₘᵢ, Pₘₐ), colormap = Reverse(:bamako), algorithm = :mip)
       Colorbar(fig[2, 4], v3.plot, label = "P (mmol N / m³)")
       v4 = volume(fig[3, 3], xnodes(grid, Center(), Center(), Center())[1:grid.Nx], ynodes(grid, Center(), Center(), Center())[1:grid.Ny], znodes(grid, Center(), Center(), Center())[1:grid.Nz], POC_plt, colorrange = (POCₘᵢ, POCₘₐ), colormap = Reverse(:bamako), algorithm = :mip)
       Colorbar(fig[3, 4], v4.plot, label = "POC (mmol C / m³)")
       record(fig, "ovs_fig.mp4", 1:length(u.times); framerate = 10) do i
           n[] = i
           msg = string("Plotting frame ", i, " of ", length(u.times))
           print(msg * " \r")
           l.text = "t=$(prettytime(u.times[i]))"
       end
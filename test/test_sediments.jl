using OceanBioME, Oceananigans, Test
using OceanBioME.Sediments: SimpleMultiG
using Oceananigans.Units

using CairoMakie

function test_flat_sediment(architecture)
    grid = RectilinearGrid(architecture; size=(3, 3, 3), extent=(10, 10, 10))

    sediment_model = SimpleMultiG(grid)

    biogeochemistry = LOBSTER(;grid, carbonates = true, oxygen = true, variable_redfield = true, open_bottom = true, sediment_model)

    model = NonhydrostaticModel(;grid, biogeochemistry, 
                                 boundary_conditions = (DIC = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂)), 
                                                        O₂ = FieldBoundaryConditions(top = GasExchange(; gas = :O₂))),
                                 tracers = (:T, :S),
                                 closure = ScalarDiffusivity(ν = 10 ^ -3, κ = 10 ^ -3))

    set!(model.biogeochemistry.sediment_model.fields.N_fast, 30.0)
    set!(model.biogeochemistry.sediment_model.fields.N_slow, 30.0)

    set!(model.biogeochemistry.sediment_model.fields.C_fast, 30.0 * model.biogeochemistry.organic_redfield)
    set!(model.biogeochemistry.sediment_model.fields.C_slow, 30.0 * model.biogeochemistry.organic_redfield)

    set!(model, P = 0.03, Z = 0.03, NO₃ = 11.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2400.0, O₂ = 240.0, 
         sPOC = model.biogeochemistry.organic_redfield, sPON = 1, bPOC = model.biogeochemistry.organic_redfield, bPON = 1,
         T = 20, S = 35)

    simulation = Simulation(model, Δt = 700.0, stop_time = 500days)

    simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                           filename = "sediment_test_tracers.jld2",
                                                           schedule = IterationInterval(1),
                                                           overwrite_existing = true)

    simulation.output_writers[:sediment] = JLD2OutputWriter(model, model.biogeochemistry.sediment_model.fields,
                                                            indices = (:, :, 1),
                                                            filename = "sediment_test_sediment.jld2",
                                                            schedule = IterationInterval(1),
                                                            overwrite_existing = true)

    @inline progress(simulation) = @info "Time: $(prettytime(simulation.model.clock.time)), Iteration: $(simulation.model.clock.iteration)"
    
    simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

    scale_negative_tracers = ScaleNegativeTracers(tracers = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON))
    simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())
    
    plankton_redfield = model.biogeochemistry.phytoplankton_redfield
    scale_negative_carbon_tracers = ScaleNegativeTracers(tracers = (:P, :Z, :DOC, :sPOC, :bPOC, :DIC), 
                                                         scalefactors = (P = plankton_redfield, 
                                                                         Z = plankton_redfield, 
                                                                         DOC = 1, sPOC = 1, bPOC = 1, DIC = 1))
    simulation.callbacks[:neg_carbon] = Callback(scale_negative_carbon_tracers; callsite = UpdateStateCallsite())

    return simulation, model
end

function load_timeseries(simulation)
    tracer_timeseries = NamedTuple()
    for tracer in keys(model.tracers)
        if !(tracer in [:T, :S])
            @info tracer
            tracer_timeseries = merge(tracer_timeseries, NamedTuple{(tracer, )}((FieldTimeSeries(simulation.output_writers[:tracers].filepath, "$tracer")[2, 2, 1:3, 1:end], )))
        end
    end

    sediment_timeseries = NamedTuple()
    for tracer in keys(model.biogeochemistry.sediment_model.fields)
        @info tracer
        sediment_timeseries = merge(sediment_timeseries, NamedTuple{(tracer, )}((FieldTimeSeries(simulation.output_writers[:sediment].filepath, "$tracer")[2, 2, 1, 1:end], )))
    end

    return tracer_timeseries, sediment_timeseries
end

function plot_timeseries(simulation, tracer_timeseries, sediment_timeseries; times = FieldTimeSeries(simulation.output_writers[:tracers].filepath, "P").times)
    fig = Figure(resolution = (1600, 1600))

    for (idx, (name, values)) in enumerate(pairs(tracer_timeseries))
        @info "$name"
        ax = Axis(fig[idx % 5, ceil(Int, idx / 5) * 2 - 1], title = "$name", xlabel = "Time (days)", ylabel = "Depth (m)")
        hm = heatmap!(ax, times ./ days, znodes(Center, simulation.model.grid), values')
        Colorbar(fig[idx % 5, ceil(Int, idx / 5) * 2], hm, label = "$name Concentration (mmol / m³)")
    end

    for (idx, (name, values)) in enumerate(pairs(sediment_timeseries))
        @info "$name"
        ax = Axis(fig[(idx + length(tracer_timeseries)) % 5, ceil(Int, (idx + length(tracer_timeseries)) / 5) * 2 - 1], title = "$name", xlabel = "Time (days)", ylabel = "Concentration (mmol / m²)")
        hm = lines!(ax, times ./ days, values)
    end

    return fig
end
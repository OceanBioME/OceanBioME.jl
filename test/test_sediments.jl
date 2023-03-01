using OceanBioME, Oceananigans, Test
using OceanBioME.Sediments: SimpleMultiG
using Oceananigans.Units

function test_flat_sediment(architecture)
    grid = RectilinearGrid(architecture; size=(3, 3, 3), extent=(10, 10, 10))

    sediment_model = SimpleMultiG(grid)

    biogeochemistry = LOBSTER(;grid, carbonates = true, oxygen = true, variable_redfield = true, open_bottom = true, sediment_model)

    model = NonhydrostaticModel(;grid, biogeochemistry, 
                                 boundary_conditions = (DIC = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂)), 
                                                        O₂ = FieldBoundaryConditions(top = GasExchange(; gas = :O₂))),
                                 tracers = (:T, :S),
                                 closure = ScalarDiffusivity(ν = 10 ^ -2, κ = 10 ^ -2))

    set!(model.biogeochemistry.sediment_model.fields.N_fast, 30.0)
    set!(model.biogeochemistry.sediment_model.fields.N_slow, 30.0)

    set!(model.biogeochemistry.sediment_model.fields.C_fast, 30.0 * model.biogeochemistry.organic_redfield)
    set!(model.biogeochemistry.sediment_model.fields.C_slow, 30.0 * model.biogeochemistry.organic_redfield)

    set!(model, P = 0.03, Z = 0.03, NO₃ = 11.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2400.0, O₂ = 240.0, 
         sPOC = model.biogeochemistry.organic_redfield, sPON = 1, bPOC = model.biogeochemistry.organic_redfield, bPON = 1,
         T = 20, S = 35)

    simulation = Simulation(model, Δt = 300.0, stop_time = 50days)

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

    return simulation, model
end
using OceanBioME, Oceananigans, Test
using OceanBioME.Sediments: Soetaert

function test_flat_sediment(architecture)
    grid = RectilinearGrid(architecture; size=(3, 3, 3), extent=(1, 1, 10))

    sediment_model = Soetaert(grid)

    biogeochemistry = LOBSTER(;grid, carbonates = true, oxygen = true, variable_redfield = true, open_bottom = true, sediment_model)

    model = NonhydrostaticModel(;grid, biogeochemistry, 
                                 boundary_conditions = (DIC = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂)), 
                                                        O₂ = FieldBoundaryConditions(top = GasExchange(; gas = :O₂))),
                                 tracers = (:T, :S))

    set!(model.biogeochemistry.sediment_model.fields.N_fast, 30.0)
    set!(model.biogeochemistry.sediment_model.fields.N_slow, 30.0)

    set!(model.biogeochemistry.sediment_model.fields.C_fast, 30.0 * model.biogeochemistry.organic_redfield)
    set!(model.biogeochemistry.sediment_model.fields.C_slow, 30.0 * model.biogeochemistry.organic_redfield)

    set!(model, P = 0.03, Z = 0.03, NO₃ = 11.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2400.0, O₂ = 240.0, 
         sPOC = model.biogeochemistry.organic_redfield, sPON = 1, bPOC = model.biogeochemistry.organic_redfield, bPON = 1,
         T = 20, S = 35)

    time_step!(model, 1.0, euler = true)
    
    return model
end
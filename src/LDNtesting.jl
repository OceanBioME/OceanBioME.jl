using OceanBioME, Test, CUDA, Oceananigans, JLD2, Documenter

using OceanBioME.Sediments: SimpleMultiG, InstantRemineralisation
using Oceananigans.Units

using OceanBioME.Sediments: sediment_tracers, sediment_fields
using Oceananigans: Field
using Oceananigans.Fields: TracerFields

using Oceananigans.Operators: volume, Azᶠᶜᶜ

using OceanBioME.LOBSTERModel: VariableRedfieldLobster



architecture = CUDA.has_cuda() ? GPU() : CPU()

function intercept_tracer_tendencies!(model, intercepted_tendencies)
    for (name, field) in enumerate(intercepted_tendencies)
        field .= Array(interior(model.timestepper.Gⁿ[name + 3]))
    end
end

function set_defaults!(sediment::SimpleMultiG)
    set!(sediment.fields.N_fast, 0.0230)
    set!(sediment.fields.N_slow, 0.0807)

    set!(sediment.fields.C_fast, 0.5893)
    set!(sediment.fields.C_slow, 0.1677)
end 



function set_defaults!(::VariableRedfieldLobster, model)

    set!(model, Z = 0.5363, 
                NO₃ = 2.3103, NH₄ = 0.0010, 
                DIC = 2106.9, Alk = 2408.9, 
                O₂ = 258.92, 
                DOC = 5.3390, DON = 0.8115,
                sPON = 0.2299, sPOC = 1.5080,
                bPON = 0.0103, bPOC = 0.0781)
    
    

    kick = 0.05
    uᵢ(x, y, z) = kick * randn()
    vᵢ(x, y, z) = kick * randn()
    wᵢ(x, y, z) = kick * randn()
    bᵢ(x, y, z) = kick * randn() + 1
    Pᵢ(x, y, z) = (1000-z)/1500
    #Pᵢ(x, y, z) = 0.4686 + exp(-((z - 500) / 50)^2)

    set!(model, u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ, P = Pᵢ)
    @info "flag"

end


total_nitrogen(sed::SimpleMultiG) = sum(sed.fields.N_fast) + 
                                    sum(sed.fields.N_slow) + 
                                    sum(sed.fields.N_ref)


total_nitrogen(::VariableRedfieldLobster, model) = sum(model.tracers.NO₃) +
                                                     sum(model.tracers.NH₄) +
                                                     sum(model.tracers.P) +
                                                     sum(model.tracers.Z) +
                                                     sum(model.tracers.DON) +
                                                     sum(model.tracers.sPON) +
                                                     sum(model.tracers.bPON)
                                                     


function test_flat_sediment(grid, biogeochemistry, model; timestepper = :QuasiAdamsBashforth2)
    Re = 5000
    model = model(; grid, 
                    biogeochemistry,
                    closure = nothing,
                    buoyancy = Buoyancy(model=BuoyancyTracer()),
                    tracers = :b)
    
                    @info "flag1"
    set_defaults!(model.biogeochemistry.sediment)

    set_defaults!(biogeochemistry.underlying_biogeochemistry, model)

    simulation = Simulation(model, Δt = 50, stop_time = 1day)
                    @info "flag2"
    intercepted_tendencies = Tuple(Array(interior(field)) for field in values(TracerFields(keys(model.tracers), grid)))

    simulation.callbacks[:intercept_tendencies] = Callback(intercept_tracer_tendencies!; callsite = TendencyCallsite(), parameters = intercepted_tendencies)

    simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "LDNtesting.jld2",
                                                       schedule = TimeInterval(24minute),
                                                       overwrite_existing = true)
@info "flag3"

    run!(simulation)
    return nothing
end

bottom_height(x, y) = -1500 + 1000 * exp(- (x^2 + y^2) / 250) # a perfect hill

display_name(::LOBSTER) = "LOBSTER"
display_name(::NutrientPhytoplanktonZooplanktonDetritus) = "NPZD"
display_name(::SimpleMultiG) = "Multi-G"
display_name(::InstantRemineralisation) = "Instant remineralisation"
display_name(::RectilinearGrid) = "Rectilinear grid"
display_name(::LatitudeLongitudeGrid) = "Latitude longitude grid"
display_name(::ImmersedBoundaryGrid) = "Immersed boundary grid"


    grid = #RectilinearGrid(architecture; size=(3, 3, 50), extent=(10, 10, 500))
    ImmersedBoundaryGrid(
        LatitudeLongitudeGrid(architecture; size = (10, 10, 50), latitude = (-10, 10), longitude = (-10, 10), z = (-1000, 0)), 
        GridFittedBottom(bottom_height))

    timestepper = :QuasiAdamsBashforth2
    sediment_model = SimpleMultiG(; grid)
    model =  HydrostaticFreeSurfaceModel
    biogeochemistry = LOBSTER(; grid,
                                carbonates = ifelse(isa(sediment_model, SimpleMultiG), true, false), 
                                oxygen = ifelse(isa(sediment_model, SimpleMultiG), true, false), 
                                variable_redfield = ifelse(isa(sediment_model, SimpleMultiG), true, false), 
                                sediment_model)
    
    @info "Running sediment on $(typeof(architecture)) with $timestepper and $(display_name(sediment_model)) on $(display_name(biogeochemistry.underlying_biogeochemistry)) with $(display_name(grid))"
    test_flat_sediment(grid, biogeochemistry, model; timestepper)                              
                
        
    





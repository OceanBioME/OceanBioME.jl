include("dependencies_for_runtests.jl")

using OceanBioME.Models: InstantRemineralisation, SimpleMultiG
using OceanBioME.Sediments: BiogeochemicalSediment

display_name(::LOBSTER) = "LOBSTER"
display_name(::NutrientPhytoplanktonZooplanktonDetritus) = "NPZD"
display_name(::BiogeochemicalSediment{<:SimpleMultiG}) = "Multi-G"
display_name(::BiogeochemicalSediment{<:InstantRemineralisation}) = "Instant remineralisation"
display_name(::RectilinearGrid) = "Rectilinear grid"
display_name(::LatitudeLongitudeGrid) = "Latitude longitude grid"
display_name(::ImmersedBoundaryGrid) = "Immersed boundary grid"

function display_name(architecture, grid, sediment_model, biogeochemistry, timestepper)
    arch_name = typeof(architecture)
    sediment_name = display_name(sediment_model)
    bgc_name = display_name(biogeochemistry.underlying_biogeochemistry)
    grid_name = display_name(grid)

    @info "Testing sediment on $arch_name with $timestepper and $sediment_name on $bgc_name with $grid_name"

    return "$architecture, $timestepper, $sediment_name, $bgc_name, $grid_name"
end

set_sinkers!(::NPZD, model) = set!(model, D = 1)
set_sinkers!(::LOBSTER, model) = set!(model, sPOM = 1, bPOM = 1)
set_sinkers!(::LOBSTER{<:Any, Val{(true, true, true)}}, model) = 
    set!(model, sPON = 1, bPON = 1, sPOC = 6.56, bPOC = 6.56)

sum_of_volume_integrals(biogeochemistry, tracers) = sum(map(f -> Field(Integral(f)), values(tracers)))
sum_of_volume_integrals(::LOBSTER{<:Any, Val{(true, true, true)}}, tracers) = 
    sum([Field(Integral(f)) for (n, f) in pairs(tracers) if n in (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)])

sum_of_area_integrals(sediment, fields) = sum(map(f -> Field(Integral(f, dims = (1, 2))), values(fields)))
sum_of_area_integrals(::SimpleMultiG{Nothing}, fields) = 
    sum([Field(Integral(f, dims = (1, 2))) for (n, f) in pairs(fields) if n in (:Nf, :Ns, :Nr)])

function test_sediment(grid, biogeochemistry, model_name, advection = WENO(order = 3, bounds = (0, 1)))
    method = quote
        return $(model_name)(; grid = $(grid), 
                               biogeochemistry = $(biogeochemistry), 
                               buoyancy = nothing, 
                               tracers = (), 
                               $(ifelse(model_name == NonhydrostaticModel, :advection, :tracer_advection)) = $advection)
    end

    model = eval(method)

    sediment_model = biogeochemistry.sediment

    set_sinkers!(biogeochemistry.underlying_biogeochemistry, model)

    if isa(biogeochemistry.sediment.biogeochemistry, SimpleMultiG)
        set!(model, NO₃ = 10, NH₄ = 1, O₂ = 1000)
    end

    tracer_nitrogen = sum_of_volume_integrals(biogeochemistry.underlying_biogeochemistry, model.tracers)

    sediment_nitrogen = sum_of_area_integrals(biogeochemistry.sediment.biogeochemistry, sediment_model.fields)

    total_nitrogen = Field(tracer_nitrogen + sediment_nitrogen)

    compute!(total_nitrogen)

    initial_total_nitrogen = CUDA.@allowscalar total_nitrogen[1, 1, 1]

    for _ in 1:100
        time_step!(model, 1)
    end

    compute!(total_nitrogen)

    final_total_nitrogen = CUDA.@allowscalar total_nitrogen[1, 1, 1]

     # simple multi-G is only good to this precision, IR is fine to default
    @test isapprox(initial_total_nitrogen, final_total_nitrogen,rtol = 0.2e-6)

    @test all(interior(sediment_nitrogen) .!= 0)

    return model
end

bottom_height(x, y) = -1000 + 500 * exp(- (x^2 + y^2) / 250) # a perfect hill

rectilinear_grid = RectilinearGrid(
    architecture;
    size = (3, 3, 50),
    extent = (10, 10, 500)
)

latlon_grid = LatitudeLongitudeGrid(
    architecture;
    size = (3, 3, 16),
    latitude = (0, 10),
    longitude = (0, 10),
    z = (-500, 0)
)

immersed_latlon_grid = ImmersedBoundaryGrid(
    latlon_grid,
    GridFittedBottom(bottom_height)
)

grids = (rectilinear_grid, latlon_grid, immersed_latlon_grid)
sediment_timesteppers = (:QuasiAdamsBashforth2, :RungeKutta3)
models = (NonhydrostaticModel, )#HydrostaticFreeSurfaceModel) # I don't think we need to test on both models anymore

@testset "Sediment integration" begin
    for grid in grids, timestepper in sediment_timesteppers
        npzd_ir = NutrientPhytoplanktonZooplanktonDetritus(; 
            grid, 
            sediment_model = InstantRemineralisationSediment(grid; timestepper)
        )

        lobster_ir = LOBSTER(;
            grid,
            sediment_model = InstantRemineralisationSediment( 
                grid;
                sinking_tracers = (:sPOM, :bPOM),
                remineralisation_reciever = :NH₄,
                timestepper
            )
        )

        simple_lobster_multi_g = LOBSTER(;
            grid,
            sediment_model = SimpleMultiGSediment(grid),
            oxygen = true
        )

        full_lobster_multi_g = LOBSTER(;
            grid,
            sediment_model = SimpleMultiGSediment(
                grid;
                sinking_nitrogen = (:sPON, :bPON),
                sinking_carbon = (:sPOC, :bPOC)
            ),
            oxygen = true,
            carbonates = true,
            variable_redfield = true
        )
        
        bgcs = [npzd_ir, lobster_ir, simple_lobster_multi_g, full_lobster_multi_g]

        for model in models, biogeochemistry in bgcs
            nonhydrostatic = (model == NonhydrostaticModel)

            grid_is_immersed = grid isa ImmersedBoundaryGrid
            grid_is_latlon = grid isa LatitudeLongitudeGrid

            # Skip incompatible combinations
            if nonhydrostatic && (grid_is_immersed || grid_is_latlon)
                continue
            end

            test_name = display_name(architecture, grid, biogeochemistry.sediment, biogeochemistry, timestepper)

            @testset "$(test_name)" begin
                test_sediment(grid, biogeochemistry, model)
            end
        end
    end
end

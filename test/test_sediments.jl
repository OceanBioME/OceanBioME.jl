include("dependencies_for_runtests.jl")

using OceanBioME.Models: InstantRemineralisation
using OceanBioME.Sediments: BiogeochemicalSediment

display_name(::LOBSTER) = "LOBSTER"
display_name(::NutrientPhytoplanktonZooplanktonDetritus) = "NPZD"
#display_name(::SimpleMultiG) = "Multi-G"
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

sum_of_integrals(tracers) = sum(map(f -> Field(Integral(f)), values(tracers)))

function test_sediment(grid, biogeochemistry, model_name, advection = WENO(order = 5, bounds = (0, 1)))
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

    tracer_nitrogen = sum_of_integrals(model.tracers)

    sediment_nitrogen = Field(Integral(sediment_model.fields.storage, dims = (1, 2)))

    total_nitrogen = Field(tracer_nitrogen + sediment_nitrogen)

    compute!(total_nitrogen)

    initial_total_nitrogen = CUDA.@allowscalar total_nitrogen[1, 1, 1]

    for _ in 1:100
        time_step!(model, 1)
    end

    compute!(total_nitrogen)

    final_total_nitrogen = CUDA.@allowscalar total_nitrogen[1, 1, 1]

    @test initial_total_nitrogen ≈ final_total_nitrogen

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

grids = (rectilinear_grid,)# latlon_grid, underlying_latlon_grid)
sediment_timesteppers = (:QuasiAdamsBashforth2, )#:RungeKutta3)
models = (NonhydrostaticModel,)# HydrostaticFreeSurfaceModel) # I don't think we need to test on both models anymore

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
        
        bgcs = [npzd_ir, lobster_ir]

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

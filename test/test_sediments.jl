include("dependencies_for_runtests.jl")

using OceanBioME.Models.SedimentModels: InstantRemineralisation, SimpleMultiG
using OceanBioME.Sediments: BiogeochemicalSediment

display_name(::NutrientsPlanktonDetritus) = "NutrientsPlanktonDetritus"
display_name(::BiogeochemicalSediment{<:Any, <:SimpleMultiG}) = "Multi-G"
display_name(::BiogeochemicalSediment{<:Any, <:InstantRemineralisation}) = "Instant remineralisation"
display_name(::RectilinearGrid) = "Rectilinear grid"
display_name(::LatitudeLongitudeGrid) = "Latitude longitude grid"
display_name(::ImmersedBoundaryGrid) = "Immersed boundary grid"

function display_name(architecture, grid, sediment_model, biogeochemistry, model_name)
    arch_name = typeof(architecture)
    sediment_name = display_name(sediment_model)
    bgc_name = display_name(biogeochemistry.underlying_biogeochemistry)
    grid_name = display_name(grid)

    @info "Testing sediment on $arch_name with $model_name and $sediment_name on $bgc_name with $grid_name"

    return "$architecture, $model_name, $sediment_name, $bgc_name, $grid_name"
end

set_sinkers!(::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:Detritus}, model) = set!(model, D = 1)
set_sinkers!(::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:DissolvedParticulate{1, 2}}, model) = set!(model, sPOM = 1, bPOM = 1)
set_sinkers!(::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonNitrogenDissolvedParticulate}, model) =
    set!(model, sPON = 1, bPON = 1, sPOC = 6.56, bPOC = 6.56)

sum_of_volume_integrals(biogeochemistry, tracers) = sum(map(f -> Field(Integral(f)), values(tracers)))
sum_of_volume_integrals(::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:CarbonNitrogenDissolvedParticulate}, tracers) =
    sum([Field(Integral(f)) for (n, f) in pairs(tracers) if n in (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)])
# `O₂` carries no nitrogen and isn't conserved against the sediment inventory, so it has to be
# excluded from what this test calls "total nitrogen"
sum_of_volume_integrals(::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:DissolvedParticulate, <:Any, <:Oxygen}, tracers) =
    sum([Field(Integral(f)) for (n, f) in pairs(tracers) if n != :O₂])

sum_of_area_integrals(sediment, fields) = sum(map(f -> Field(Integral(f, dims = (1, 2))), values(fields)))
sum_of_area_integrals(::SimpleMultiG{Nothing}, fields) =
    sum([Field(Integral(f, dims = (1, 2))) for (n, f) in pairs(fields) if n in (:Nf, :Ns, :Nr)])

function test_sediment(grid, biogeochemistry, model_name, advection = WENO(order = 3, bounds = (0, 1)))
    method = quote
        return $(model_name)($(grid);
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

    @test isapprox(initial_total_nitrogen, final_total_nitrogen, rtol = 1e-8)#0.2e-6)

    @test CUDA.@allowscalar all(interior(sediment_nitrogen) .!= 0)

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
models = (NonhydrostaticModel, HydrostaticFreeSurfaceModel) # exercises both `substep_sediment!` methods (RK3 and AB2)

@testset "Sediment integration" begin
    for grid in grids
        npzd_ir = NPZD(
            grid;
            sediment = InstantRemineralisationSediment(grid)
        )

        lobster_ir = LOBSTER(
            grid;
            sediment = InstantRemineralisationSediment(
                grid;
                sinking_tracers = (:sPOM, :bPOM),
                remineralisation_reciever = :NH₄
            )
        )

        simple_lobster_multi_g = LOBSTER(
            grid;
            sediment = SimpleMultiGSediment(grid),
            oxygen = Oxygen()
        )

        full_lobster_multi_g = LOBSTER(
            grid;
            detritus = CarbonNitrogenDissolvedParticulate(grid; open_bottom = true),
            sediment = SimpleMultiGSediment(
                grid;
                sinking_nitrogen = (:sPON, :bPON),
                sinking_carbon = (:sPOC, :bPOC)
            ),
            oxygen = Oxygen(),
            inorganic_carbon = CarbonateSystem()
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

            test_name = display_name(architecture, grid, biogeochemistry.sediment, biogeochemistry, model)

            @testset "$(test_name)" begin
                # `InstantRemineralisation` defines its tendency method for the chosen
                # `remineralisation_reciever` via `eval` at construction time, so a top-level
                # loop that both builds the biogeochemistry and calls `test_sediment` in the
                # same compiled thunk needs `invokelatest` to see it
                Base.invokelatest(test_sediment, grid, biogeochemistry, model)
            end
        end
    end
end

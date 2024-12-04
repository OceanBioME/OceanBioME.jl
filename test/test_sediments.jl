include("dependencies_for_runtests.jl")

using Oceananigans.Units

using Oceananigans: Field, TendencyCallsite
using Oceananigans.Fields: TracerFields
using Oceananigans.Operators: volume, Azᶠᶜᶜ

using OceanBioME.Sediments: SimpleMultiG, InstantRemineralisation, sediment_tracers, sediment_fields
using OceanBioME.Models.LOBSTERModel: VariableRedfieldLobster

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

set_defaults!(::InstantRemineralisation) = nothing

set_defaults!(::LOBSTER, model) =
    set!(model, P = 0.4686, Z = 0.5363,
                NO₃ = 2.3103, NH₄ = 0.0010,
                DOM = 0.8115,
                sPOM = 0.2299, bPOM = 0.0103)


set_defaults!(::VariableRedfieldLobster, model) =
    set!(model, P = 0.4686, Z = 0.5363,
                NO₃ = 2.3103, NH₄ = 0.0010,
                DIC = 2106.9, Alk = 2408.9,
                O₂ = 258.92,
                DOC = 5.3390, DON = 0.8115,
                sPON = 0.2299, sPOC = 1.5080,
                bPON = 0.0103, bPOC = 0.0781)


set_defaults!(::NutrientPhytoplanktonZooplanktonDetritus, model) =  set!(model, N = 2.3, P = 0.4, Z = 0.5, D = 0.2)

total_nitrogen(sed::SimpleMultiG) = sum(sed.fields.N_fast) +
                                    sum(sed.fields.N_slow) +
                                    sum(sed.fields.N_ref)

total_nitrogen(sed::InstantRemineralisation) = sum(sed.fields.N_storage)

total_nitrogen(::LOBSTER, model) = sum(model.tracers.NO₃) +
                                   sum(model.tracers.NH₄) +
                                   sum(model.tracers.P) +
                                   sum(model.tracers.Z) +
                                   sum(model.tracers.DOM) +
                                   sum(model.tracers.sPOM) +
                                   sum(model.tracers.bPOM)

total_nitrogen(::VariableRedfieldLobster, model) = sum(model.tracers.NO₃) +
                                                     sum(model.tracers.NH₄) +
                                                     sum(model.tracers.P) +
                                                     sum(model.tracers.Z) +
                                                     sum(model.tracers.DON) +
                                                     sum(model.tracers.sPON) +
                                                     sum(model.tracers.bPON)

total_nitrogen(::NutrientPhytoplanktonZooplanktonDetritus, model) = sum(model.tracers.N) +
                                                                    sum(model.tracers.P) +
                                                                    sum(model.tracers.Z) +
                                                                    sum(model.tracers.D)

function test_flat_sediment(grid, biogeochemistry, model; timestepper = :QuasiAdamsBashforth2)
    model = isa(model, NonhydrostaticModel) ? model(; grid,
                                                      biogeochemistry,
                                                      closure = nothing,
                                                      timestepper,
                                                      buoyancy = nothing) :
                                              model(; grid,
                                                      biogeochemistry,
                                                      closure = nothing,
                                                      buoyancy = nothing,
                                                      tracers = nothing)

    set_defaults!(model.biogeochemistry.sediment)

    set_defaults!(biogeochemistry.underlying_biogeochemistry, model)

    simulation = Simulation(model, Δt = 50, stop_time = 1day)

    intercepted_tendencies = Tuple(Array(interior(field)) for field in values(TracerFields(keys(model.tracers), grid)))

    simulation.callbacks[:intercept_tendencies] = Callback(intercept_tracer_tendencies!; callsite = TendencyCallsite(), parameters = intercepted_tendencies)

    N₀ = CUDA.@allowscalar total_nitrogen(biogeochemistry.underlying_biogeochemistry, model) * volume(1, 1, 1, grid, Center(), Center(), Center()) + total_nitrogen(biogeochemistry.sediment) * Azᶠᶜᶜ(1, 1, 1, grid)

    run!(simulation)

    # the model is changing the tracer tendencies
    @test any([any(intercepted_tendencies[idx] .!= Array(interior(model.timestepper.Gⁿ[tracer]))) for (idx, tracer) in enumerate(keys(model.tracers))])

    # the sediment tendencies are being updated
    @test all([any(Array(interior(tend)) .!= 0.0) for tend in model.biogeochemistry.sediment.tendencies.Gⁿ])
    @test all([any(Array(interior(tend)) .!= 0.0) for tend in model.biogeochemistry.sediment.tendencies.G⁻])

    # the sediment values are being integrated
    initial_values = (N_fast = 0.0230, N_slow = 0.0807, C_fast = 0.5893, C_slow = 0.1677, N_ref = 0.0, C_ref = 0.0, N_storage = 0.0)
    @test all([any(Array(interior(field)) .!= initial_values[name]) for (name, field) in pairs(model.biogeochemistry.sediment.fields)])

    N₁ = CUDA.@allowscalar total_nitrogen(biogeochemistry.underlying_biogeochemistry, model) * volume(1, 1, 1, grid, Center(), Center(), Center()) + total_nitrogen(biogeochemistry.sediment) * Azᶠᶜᶜ(1, 1, 1, grid)

    # conservations
    rtol = ifelse(isa(architecture, CPU), max(√eps(N₀), √eps(N₁)), 5e-7)
    @test isapprox(N₀, N₁; rtol)

    return nothing
end

display_name(::LOBSTER) = "LOBSTER"
display_name(::NutrientPhytoplanktonZooplanktonDetritus) = "NPZD"
display_name(::SimpleMultiG) = "Multi-G"
display_name(::InstantRemineralisation) = "Instant remineralisation"
display_name(::RectilinearGrid) = "Rectilinear grid"
display_name(::LatitudeLongitudeGrid) = "Latitude longitude grid"
display_name(::ImmersedBoundaryGrid) = "Immersed boundary grid"


bottom_height(x, y) = -1000 + 500 * exp(- (x^2 + y^2) / 250) # a perfect hill

@testset "Sediment integration" begin
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

    underlying_latlon_grid = LatitudeLongitudeGrid(
        architecture;
        size = (3, 3, 16),
        latitude = (0, 10),
        longitude = (0, 10),
        z = (-500, 0)
    )

    immersed_latlon_grid = ImmersedBoundaryGrid(
        underlying_latlon_grid,
        GridFittedBottom(bottom_height)
    )

    grids = (rectilinear_grid, latlon_grid, underlying_latlon_grid)
    timesteppers = (:QuasiAdamsBashforth2, :RungeKutta3)
    sediment_models = (InstantRemineralisation(; grid), SimpleMultiG(; grid))
    models = (NonhydrostaticModel, HydrostaticFreeSurfaceModel)

    for grid in grids, timestepper in timesteppers, sediment_model in sediment_models, model in models
        npzd_bgc_model = NutrientPhytoplanktonZooplanktonDetritus(; grid, sediment_model)

        using_simple_multi_g = sediment_model isa SimpleMultiG

        lobster_bgc_model = LOBSTER(;
            grid,
            sediment_model,
            carbonates = using_simple_multi_g,
            oxygen = using_simple_multi_g,
            variable_redfield = using_simple_multi_g,
        )

        bgc_models = [npzd_bgc_model, lobster_bgc_model]

        for biogeochemistry in bgc_models
            nonhydrostatic = (model == NonhydrostaticModel)
            hydrostatic = (model == HydrostaticFreeSurfaceModel)

            grid_is_immersed = grid isa ImmersedBoundaryGrid
            grid_is_latlon = grid isa LatitudeLongitudeGrid

            rk3 = (timestepper == :RungeKutta3)
            simple_multi_g = sediment_model isa SimpleMultiG
            bgc_is_npzd = biogeochemistry.underlying_biogeochemistry isa NutrientPhytoplanktonZooplanktonDetritus

            # Skip incompatible combinations
            skip1 = nonhydrostatic && (grid_is_immersed || grid_is_latlon)
            skip2 = hydrostatic && rk3
            skip3 = simple_multi_g && bgc_is_npzd

            if skip1 || skip2 || skip3
                continue
            end

            arch_name = typeof(architecture)
            sediment_name = display_name(sediment_model)
            bgc_name = display_name(biogeochemistry.underlying_biogeochemistry)
            grid_name = display_name(grid)

            @info "Testing sediment on $arch_name with $timestepper and $sediment_name on $bgc_name with $grid_name"

            @testset "$architecture, $timestepper, $sediment_name, $bgc_name, $grid_name" begin
                @test_skip test_flat_sediment(grid, biogeochemistry, model; timestepper)
            end
        end
    end
end

include("dependencies_for_runtests.jl")

using CUDA: @allowscalar

using OceanBioME: TwoBandPhotosyntheticallyActiveRadiation, LOBSTER, NPZD

using Oceananigans.Architectures: on_architecture
using Oceananigans.Biogeochemistry: update_biogeochemical_state!, required_biogeochemical_tracers, biogeochemical_auxiliary_fields

Pᵢ(x,y,z) = 2.5 + z

function test_two_band(grid, model_type, surface_PAR, discrete_form)
    biogeochemistry = NPZD(grid; 
                           light_attenuation = 
                               TwoBandPhotosyntheticallyActiveRadiation(; grid,
                                                                       surface_PAR, 
                                                                       discrete_form))

    model = model_type(grid;
                       biogeochemistry,
                       tracers = unique((required_biogeochemical_tracers(biogeochemistry)..., :T, :S))) # because hydrostatic free surface will request T and S and some BGC models will too

    set!(model, P = Pᵢ)

    PAR_model = model.biogeochemistry.light_attenuation

    kʳ = PAR_model.water_red_attenuation
    kᵇ = PAR_model.water_blue_attenuation
    χʳ = PAR_model.chlorophyll_red_attenuation
    χᵇ = PAR_model.chlorophyll_blue_attenuation
    eʳ = PAR_model.chlorophyll_red_exponent
    eᵇ = PAR_model.chlorophyll_blue_exponent
    r = PAR_model.pigment_ratio
    Rᶜₚ = PAR_model.phytoplankton_chlorophyll_ratio

    ∫Chlʳ = [(2.0 * Rᶜₚ / r) ^ eʳ * 0.5]
    ∫Chlᵇ = [(2.0 * Rᶜₚ / r) ^ eᵇ * 0.5]

    push!(∫Chlʳ, ∫Chlʳ[1] + (2.0 * Rᶜₚ / r) ^ eʳ * 0.5 + (1.0 * Rᶜₚ / r) ^ eʳ * 0.5)
    push!(∫Chlᵇ, ∫Chlᵇ[1] + (2.0 * Rᶜₚ / r) ^ eᵇ * 0.5 + (1.0 * Rᶜₚ / r) ^ eᵇ * 0.5)

    expected_PAR = 100.0 .* [exp(- 0.5 * kʳ - ∫Chlʳ[1] * χʳ) + exp(- 0.5 * kᵇ - ∫Chlᵇ[1] * χᵇ),
                             exp(- 1.5 * kʳ - ∫Chlʳ[2] * χʳ) + exp(- 1.5 * kᵇ - ∫Chlᵇ[2] * χᵇ)] ./ 2

    results_PAR = Array(interior(biogeochemical_auxiliary_fields(biogeochemistry).PAR))[1, 1, 1:2]

    @test all(results_PAR .≈ reverse(expected_PAR))

    return nothing
end

function test_multi_band(grid, model_type, surface_PAR, discrete_form)
    light_attenuation = MultiBandPhotosyntheticallyActiveRadiation(; grid,
                                                                     bands = ((1, 2), ),
                                                                     base_bands = [1, 2],
                                                                     base_water_attenuation_coefficient = [0.01, 0.01],
                                                                     base_chlorophyll_exponent = [2, 2],
                                                                     base_chlorophyll_attenuation_coefficient = [0.1, 0.1],
                                                                     surface_PAR, discrete_form)

    biogeochemistry = NPZD(grid; light_attenuation)

    model = model_type(grid;
                       biogeochemistry,
                       buoyancy = nothing,
                       tracers = nothing)

    set!(model, P = 2/1.31) # this will cause tests to fail for models with different chlorophyll ratios

    expected_PAR = on_architecture(CPU(), 100 .* exp.(znodes(grid, Center()) * (0.01 + 0.1 * 2 ^ 2)))

    @test (@allowscalar all(interior(on_architecture(CPU(), light_attenuation.fields[1]), 1, 1, :) .≈ expected_PAR))

    light_attenuation = MultiBandPhotosyntheticallyActiveRadiation(; grid,
                                                                     bands = ((1, 2), (8, 9)),
                                                                     base_bands = [1, 2, 8, 9],
                                                                     base_water_attenuation_coefficient = [0.01, 0.01, 0.02, 0.02],
                                                                     base_chlorophyll_exponent = [2, 2, 1.5, 1.5],
                                                                     base_chlorophyll_attenuation_coefficient = [0.1, 0.1, 0.2, 0.2],
                                                                     surface_PAR,
                                                                     discrete_form)

    biogeochemistry = NPZD(grid; light_attenuation)

    model = model_type(grid;
                       biogeochemistry,
                       buoyancy = nothing,
                       tracers = nothing)

    set!(model, P = 2/1.31) # this will cause tests to fail for models with different chlorophyll ratios (e.g. PISCES)

    expected_PAR1 = on_architecture(CPU(), 100 .* exp.(znodes(grid, Center()) * (0.01 + 0.1 * 2 ^ 2)) / 2)
    expected_PAR2 = on_architecture(CPU(), 100 .* exp.(znodes(grid, Center()) * (0.02 + 0.2 * 2 ^ 1.5)) / 2)

    PAR, PAR₁, PAR₂ = map(v-> on_architecture(CPU(), v), values(biogeochemical_auxiliary_fields(light_attenuation)))

    # not sure why I've had to reduce the tolerances here
    @test all(isapprox.(interior(PAR₁, 1, 1, :), expected_PAR1, atol=1e-4))
    @test all(isapprox.(interior(PAR₂, 1, 1, :), expected_PAR2, atol=1e-4))
    @test  all(isapprox.(PAR[1, 1, 1:grid.Nz], expected_PAR1 .+ expected_PAR2, atol=1e-3)) # binary operation so we can't `interior` it

    # check all the models work as expected
    @test isnothing(time_step!(model, 1))

    return nothing
end

@inline discrete_surface_PAR(i, j, grid, clock, fields) = 100
@inline continuous_surface_PAR(x, y, t) = 100
field_surface_PAR = Oceananigans.Fields.ConstantField(100)

@testset "Light attenuaiton model" begin
    for model in (NonhydrostaticModel, HydrostaticFreeSurfaceModel),
        grid in (RectilinearGrid(architecture; size = (2, 2, 2), extent = (2, 2, 2)),
                 LatitudeLongitudeGrid(architecture; size = (5, 5, 2), longitude = (-180, 180), latitude = (-85, 85), z = (-2, 0)))

        if !((model == NonhydrostaticModel) && ((grid isa LatitudeLongitudeGrid) | (grid isa OrthogonalSphericalShellGrid)))
            @info "Testing light with in $model on $grid..."
            test_two_band(grid, model, field_surface_PAR, false)
            test_multi_band(grid, model, field_surface_PAR, false)
        end
    end

    grid = RectilinearGrid(architecture; size = (2, 2, 2), extent = (2, 2, 2))

    for surface_PAR in (discrete_surface_PAR, continuous_surface_PAR, field_surface_PAR)
        discrete_form = surface_PAR == discrete_surface_PAR

        test_two_band(grid, NonhydrostaticModel, surface_PAR, discrete_form)
        test_multi_band(grid, NonhydrostaticModel, surface_PAR, discrete_form)
    end
end

@testset "Float32 TwoBandPhotosyntheticallyActiveRadiation" begin
    grid = RectilinearGrid(architecture, Float32; size=(3, 3, 10), extent=(10, 10, 200))
    par = TwoBandPhotosyntheticallyActiveRadiation(; grid)

    @test par.water_red_attenuation isa Float32
    @test par.water_blue_attenuation isa Float32
    @test par.chlorophyll_red_attenuation isa Float32
    @test par.chlorophyll_blue_attenuation isa Float32
    @test par.chlorophyll_red_exponent isa Float32
    @test par.chlorophyll_blue_exponent isa Float32
    @test par.pigment_ratio isa Float32
    @test par.phytoplankton_chlorophyll_ratio isa Float32
end

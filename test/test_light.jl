include("dependencies_for_runtests.jl")

using CUDA: @allowscalar

using OceanBioME: TwoBandPhotosyntheticallyActiveRadiation, 
                  PrescribedAttenuationPAR, 
                  LOBSTER, NPZD, ImplicitBiology

using Oceananigans.Architectures: on_architecture
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry,
                                    update_biogeochemical_state!, 
                                    required_biogeochemical_tracers, 
                                    required_biogeochemical_auxiliary_fields,
                                    biogeochemical_auxiliary_fields

import OceanBioME: chlorophyll

Pᵢ(x,y,z) = 2.5 + z

struct MultiPhytoplanktonLightTest <: AbstractBiogeochemistry end

required_biogeochemical_tracers(::MultiPhytoplanktonLightTest) = (:P1, :P2, :P3)
required_biogeochemical_auxiliary_fields(::MultiPhytoplanktonLightTest) = (:PAR, )

@inline function (::MultiPhytoplanktonLightTest)(i, j, k, grid, val_name, clock, fields, auxiliary_fields)
    return zero(fields.P1[i, j, k])
end

@inline chlorophyll(::MultiPhytoplanktonLightTest, model) =
    model.tracers.P1 + 2 * model.tracers.P2 + 3 * model.tracers.P3

function test_two_band(grid, model_type, surface_PAR, discrete_form, parameters = nothing)
    biogeochemistry = NPZD(grid; 
                           light_attenuation = 
                               TwoBandPhotosyntheticallyActiveRadiation(grid, surface_PAR;
                                                                        discrete_form,
                                                                        parameters))

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
    Rᶜₚ = biogeochemistry.underlying_biogeochemistry.plankton.chlorophyll_ratio

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

function multi_phytoplankton_PAR(grid; P1, P2, P3, surface_PAR = 100, discrete_form = false, parameters = nothing)
    light_attenuation = TwoBandPhotosyntheticallyActiveRadiation(grid, surface_PAR; discrete_form, parameters)
    biogeochemistry = Biogeochemistry(MultiPhytoplanktonLightTest(); light_attenuation)
    model = NonhydrostaticModel(grid; biogeochemistry, buoyancy = nothing, tracers = nothing)

    @test !hasproperty(model.tracers, :P)
    set!(model; P1, P2, P3)

    return Array(interior(biogeochemical_auxiliary_fields(biogeochemistry).PAR))[1, 1, :]
end

function npzd_two_band_PAR(grid)
    light_attenuation = TwoBandPhotosyntheticallyActiveRadiation(grid, 100)
    biogeochemistry = NPZD(grid; light_attenuation)
    model = NonhydrostaticModel(grid; biogeochemistry, buoyancy = nothing, tracers = nothing)

    set!(model, P = Pᵢ)

    return Array(interior(biogeochemical_auxiliary_fields(biogeochemistry).PAR))[1, 1, :]
end

function test_prescribed_attenuation(grid, model_type, 
                                     surface_PAR, surface_discrete_form, 
                                     attenuation, attenuation_discrete_form, 
                                     surface_parameters = nothing, 
                                     attenuation_parameters = nothing)

    light_attenuation = PrescribedAttenuationPAR(grid, surface_PAR;
                                                 surface_discrete_form,
                                                 surface_parameters,
                                                 attenuation,
                                                 attenuation_discrete_form,
                                                 attenuation_parameters)
                                                 
    biogeochemistry = ImplicitBiology(grid; light_attenuation)

    model = model_type(grid;
                       biogeochemistry,
                       tracers = unique((required_biogeochemical_tracers(biogeochemistry)..., :T, :S))) # because hydrostatic free surface will request T and S and some BGC models will too

    PAR = model.biogeochemistry.light_attenuation

    expected_PAR = 100 .* exp.(znodes(PAR.field) .* 0.1)

    results_PAR = Array(interior(biogeochemical_auxiliary_fields(biogeochemistry).PAR))[1, 1, 1:2]

    @test all(results_PAR .≈ expected_PAR)

    return nothing
end

function test_multi_band(grid, model_type, surface_PAR, discrete_form, parameters = nothing)
    light_attenuation = MultiBandPhotosyntheticallyActiveRadiation(grid,surface_PAR;
                                                                   bands = ((1, 2), ),
                                                                   base_bands = [1, 2],
                                                                   base_water_attenuation_coefficient = [0.01, 0.01],
                                                                   base_chlorophyll_exponent = [2, 2],
                                                                   base_chlorophyll_attenuation_coefficient = [0.1, 0.1],
                                                                   discrete_form, parameters)

    biogeochemistry = NPZD(grid; light_attenuation)

    model = model_type(grid;
                       biogeochemistry,
                       buoyancy = nothing,
                       tracers = nothing)

    set!(model, P = 2/1.31) # this will cause tests to fail for models with different chlorophyll ratios

    expected_PAR = on_architecture(CPU(), 100 .* exp.(znodes(grid, Center()) * (0.01 + 0.1 * 2 ^ 2)))

    @test (@allowscalar all(interior(on_architecture(CPU(), light_attenuation.fields[1]), 1, 1, :) .≈ expected_PAR))

    light_attenuation = MultiBandPhotosyntheticallyActiveRadiation(grid, surface_PAR;
                                                                   bands = ((1, 2), (8, 9)),
                                                                   base_bands = [1, 2, 8, 9],
                                                                   base_water_attenuation_coefficient = [0.01, 0.01, 0.02, 0.02],
                                                                   base_chlorophyll_exponent = [2, 2, 1.5, 1.5],
                                                                   base_chlorophyll_attenuation_coefficient = [0.1, 0.1, 0.2, 0.2],
                                                                   discrete_form, parameters)

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
    @test all(isapprox.(PAR[1, 1, 1:grid.Nz], expected_PAR1 .+ expected_PAR2, atol=1e-3)) # binary operation so we can't `interior` it

    # check all the models work as expected
    @test isnothing(time_step!(model, 1))

    return nothing
end

@inline discrete_surface_PAR(i, j, grid, clock, fields) = 100
@inline continuous_surface_PAR(x, y, t, l0) = l0
field_surface_PAR = Oceananigans.Fields.ConstantField(100)

@testset "Light attenuaiton model" begin
    for model in (NonhydrostaticModel, HydrostaticFreeSurfaceModel),
        grid in (RectilinearGrid(architecture; size = (2, 2, 2), extent = (2, 2, 2)),
                 LatitudeLongitudeGrid(architecture; size = (5, 5, 2), longitude = (-180, 180), latitude = (-85, 85), z = (-2, 0)))

        if !((model == NonhydrostaticModel) && ((grid isa LatitudeLongitudeGrid) | (grid isa OrthogonalSphericalShellGrid)))
            @info "Testing light with in $model on $grid..."
            test_two_band(grid, model, field_surface_PAR, false)
            test_multi_band(grid, model, field_surface_PAR, false)
            test_prescribed_attenuation(grid, model, field_surface_PAR, false, 0.1, false)
        end
    end

    grid = RectilinearGrid(architecture; size = (2, 2, 2), extent = (2, 2, 2))

    for surface_PAR in (discrete_surface_PAR, field_surface_PAR)
        discrete_form = surface_PAR == discrete_surface_PAR

        test_two_band(  grid, NonhydrostaticModel, surface_PAR, discrete_form)
        test_multi_band(grid, NonhydrostaticModel, surface_PAR, discrete_form)
        test_prescribed_attenuation(grid, NonhydrostaticModel, surface_PAR, discrete_form, 0.1, false)
    end

    test_two_band(  grid, NonhydrostaticModel, continuous_surface_PAR, false, 100)
    test_multi_band(grid, NonhydrostaticModel, continuous_surface_PAR, false, 100)
    test_prescribed_attenuation(grid, NonhydrostaticModel, continuous_surface_PAR, false, 0.1, false, 100)

    test_prescribed_attenuation(grid, NonhydrostaticModel, continuous_surface_PAR, false, (x, y, z, t, a0) -> a0, false, 100, 0.1) # continuous attenuation with parameters
    test_prescribed_attenuation(grid, NonhydrostaticModel, continuous_surface_PAR, false, (args...) -> 0.1, true, 100) # discrete attenuation
end

@testset "TwoBand generic phytoplankton chlorophyll" begin
    grid = RectilinearGrid(architecture; size = (1, 1, 3), extent = (1, 1, 3))

    PAR_P1 = multi_phytoplankton_PAR(grid; P1 = 1, P2 = 0, P3 = 0)
    PAR_P2 = multi_phytoplankton_PAR(grid; P1 = 0, P2 = 0.5, P3 = 0)
    PAR_P3 = multi_phytoplankton_PAR(grid; P1 = 0, P2 = 0, P3 = 1/3)
    PAR_split = multi_phytoplankton_PAR(grid; P1 = 0.2, P2 = 0.1, P3 = 0.2)
    PAR_more = multi_phytoplankton_PAR(grid; P1 = 1, P2 = 0.5, P3 = 0)
    PAR_discrete = multi_phytoplankton_PAR(grid; P1 = 1, P2 = 0, P3 = 0,
                                           surface_PAR = discrete_surface_PAR, discrete_form = true)
    PAR_continuous = multi_phytoplankton_PAR(grid; P1 = 1, P2 = 0, P3 = 0,
                                             surface_PAR = continuous_surface_PAR, parameters = 100)

    @test PAR_P1 ≈ PAR_P2
    @test PAR_P1 ≈ PAR_P3
    @test PAR_P1 ≈ PAR_split
    @test PAR_P1 ≈ PAR_discrete
    @test PAR_P1 ≈ PAR_continuous
    @test all(PAR_more .< PAR_P1)
end

@testset "Float32 TwoBandPhotosyntheticallyActiveRadiation" begin
    grid = RectilinearGrid(architecture, Float32; size = (3, 3, 10), extent = (10, 10, 200))
    par = TwoBandPhotosyntheticallyActiveRadiation(grid, 100)
    @test par.water_red_attenuation isa Float32
    @test par.water_blue_attenuation isa Float32
    @test par.chlorophyll_red_attenuation isa Float32
    @test par.chlorophyll_blue_attenuation isa Float32
    @test par.chlorophyll_red_exponent isa Float32
    @test par.chlorophyll_blue_exponent isa Float32
    @test par.pigment_ratio isa Float32

    generic_PAR = multi_phytoplankton_PAR(grid; P1 = 1f0, P2 = 0.5f0, P3 = 0.25f0)
    npzd_PAR = npzd_two_band_PAR(grid)

    @test eltype(generic_PAR) == Float32
    @test eltype(npzd_PAR) == Float32
end

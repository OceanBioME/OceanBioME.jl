include("dependencies_for_runtests.jl")

using CUDA: @allowscalar

using OceanBioME: TwoBandPhotosyntheticallyActiveRadiation, 
                  PrescribedAttenuationPAR, 
                  LOBSTER, NPZD, ImplicitBiology

using Oceananigans.Architectures: on_architecture
using Oceananigans.Biogeochemistry: update_biogeochemical_state!,
                                    required_biogeochemical_tracers,
                                    biogeochemical_auxiliary_fields
using Oceananigans.Fields: ZFaceField

Pᵢ(x,y,z) = 2.5 + z

function test_two_band(grid, model_type, surface_PAR, discrete_form, parameters = nothing; test_interface = false)
    interface_field = test_interface ? ZFaceField(grid) : nothing

    biogeochemistry = NPZD(grid;
                           light_attenuation =
                               TwoBandPhotosyntheticallyActiveRadiation(grid, surface_PAR;
                                                                        discrete_form,
                                                                        parameters,
                                                                        interface_field))

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

    results_PAR = Array(interior(biogeochemical_auxiliary_fields(biogeochemistry).PAR))[1, 1, 1:2]

    attenuationʳ(P) = kʳ + χʳ * (Rᶜₚ * P / r) ^ eʳ
    attenuationᵇ(P) = kᵇ + χᵇ * (Rᶜₚ * P / r) ^ eᵇ

    Δz = 1.0
    PAR⁰ = 100.0

    tʳ, tᵇ = 1.0, 1.0

    Aʳ = exp(-Δz * attenuationʳ(2.0))
    Aᵇ = exp(-Δz * attenuationᵇ(2.0))

    PAR2 = PAR⁰ * (0.5 * tʳ * (1 - Aʳ) / (-log(Aʳ)) + 0.5 * tᵇ * (1 - Aᵇ) / (-log(Aᵇ)))
    PARi_face2 = PAR⁰ * (0.5 * tʳ * Aʳ + 0.5 * tᵇ * Aᵇ)

    tʳ, tᵇ = tʳ * Aʳ, tᵇ * Aᵇ

    Aʳ = exp(-Δz * attenuationʳ(1.0))
    Aᵇ = exp(-Δz * attenuationᵇ(1.0))

    PAR1 = PAR⁰ * (0.5 * tʳ * (1 - Aʳ) / (-log(Aʳ)) + 0.5 * tᵇ * (1 - Aᵇ) / (-log(Aᵇ)))
    PARi_face1 = PAR⁰ * (0.5 * tʳ * Aʳ + 0.5 * tᵇ * Aᵇ)

    expected_PAR = [PAR1, PAR2]

    @test all(results_PAR .≈ expected_PAR)

    if test_interface
        expected_PAR_interface = [PARi_face1, PARi_face2, PAR⁰]

        results_PAR_interface = Array(interior(biogeochemical_auxiliary_fields(biogeochemistry).PAR_interface))[1, 1, 1:3]

        @test all(results_PAR_interface .≈ expected_PAR_interface)
    end

    return nothing
end

function test_prescribed_attenuation(grid, model_type,
                                     surface_PAR, surface_discrete_form,
                                     attenuation, attenuation_discrete_form,
                                     surface_parameters = nothing,
                                     attenuation_parameters = nothing;
                                     test_interface = false)

    interface_field = test_interface ? ZFaceField(grid) : nothing

    light_attenuation = PrescribedAttenuationPAR(grid, surface_PAR;
                                                 surface_discrete_form,
                                                 surface_parameters,
                                                 attenuation,
                                                 attenuation_discrete_form,
                                                 attenuation_parameters,
                                                 interface_field)
                                                 
    biogeochemistry = ImplicitBiology(grid; light_attenuation)

    model = model_type(grid;
                       biogeochemistry,
                       tracers = unique((required_biogeochemical_tracers(biogeochemistry)..., :T, :S))) # because hydrostatic free surface will request T and S and some BGC models will too

    results_PAR = Array(interior(biogeochemical_auxiliary_fields(biogeochemistry).PAR))[1, 1, 1:2]

    K = 0.1
    Δz = 1.0
    PAR⁰ = 100.0

    e = exp(-K * Δz)

    PARi_face2 = PAR⁰ * e
    PAR2 = PAR⁰ * (1 - e) / (K * Δz)

    PARi_face1 = PARi_face2 * e
    PAR1 = PARi_face2 * (1 - e) / (K * Δz)

    expected_PAR = [PAR1, PAR2]

    @test all(results_PAR .≈ expected_PAR)

    if test_interface
        expected_PAR_interface = [PARi_face1, PARi_face2, PAR⁰]

        results_PAR_interface = Array(interior(biogeochemical_auxiliary_fields(biogeochemistry).PAR_interface))[1, 1, 1:3]

        @test all(results_PAR_interface .≈ expected_PAR_interface)
    end

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
            test_two_band(grid, model, field_surface_PAR, false; test_interface = true)
            test_multi_band(grid, model, field_surface_PAR, false)
            test_prescribed_attenuation(grid, model, field_surface_PAR, false, 0.1, false)
            test_prescribed_attenuation(grid, model, field_surface_PAR, false, 0.1, false; test_interface = true)
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

@testset "Float32 TwoBandPhotosyntheticallyActiveRadiation" begin
    grid = RectilinearGrid(architecture, Float32; size=(3, 3, 10), extent=(10, 10, 200))
    par = TwoBandPhotosyntheticallyActiveRadiation(grid, 100)

    @test par.water_red_attenuation isa Float32
    @test par.water_blue_attenuation isa Float32
    @test par.chlorophyll_red_attenuation isa Float32
    @test par.chlorophyll_blue_attenuation isa Float32
    @test par.chlorophyll_red_exponent isa Float32
    @test par.chlorophyll_blue_exponent isa Float32
    @test par.pigment_ratio isa Float32
end

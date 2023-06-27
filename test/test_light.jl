using Oceananigans, Test
using OceanBioME: TwoBandPhotosyntheticallyActiveRadiation, LOBSTER
using Oceananigans.Biogeochemistry: update_biogeochemical_state!

grid = RectilinearGrid(size=(1,1,2), extent=(1,1,2))

@testset "Two band attenuation" begin

    model = NonhydrostaticModel(; grid, 
                                  biogeochemistry = LOBSTER(; grid,
                                                              light_attenuation_model = TwoBandPhotosyntheticallyActiveRadiation(; grid),
                                                              surface_phytosynthetically_active_radiation = (x, y, t) -> 100.0))
    Pᵢ(x,y,z) = 2.5 + z

    set!(model, P = Pᵢ)

    update_biogeochemical_state!(model.biogeochemistry, model)

    PAR_model = model.biogeochemistry.light_attenuation_model
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

    results_PAR = convert(Array, model.biogeochemistry.light_attenuation_model.field)[1, 1, 1:2]

    @test all(results_PAR .≈ reverse(expected_PAR))
end

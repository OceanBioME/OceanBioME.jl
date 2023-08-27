using Oceananigans, Test
using OceanBioME: TwoBandPhotosyntheticallyActiveRadiation, LOBSTER, NutrientPhytoplanktonZooplanktonDetritus
using Oceananigans.Biogeochemistry: update_biogeochemical_state!, required_biogeochemical_tracers, biogeochemical_auxiliary_fields

Pᵢ(x,y,z) = 2.5 + z

function test_two_band(grid, bgc, model_type)
    biogeochemistry = bgc(; grid,
                            light_attenuation_model = TwoBandPhotosyntheticallyActiveRadiation(; grid),
                            surface_phytosynthetically_active_radiation = (x, y, t) -> 100.0)

    model = model_type(; grid, 
                         biogeochemistry,
                         tracers = unique((required_biogeochemical_tracers(biogeochemistry)..., :T, :S))) # because hydrostatic free surface will request T and S and some BGC models will too

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

    results_PAR = convert(Array, biogeochemical_auxiliary_fields(biogeochemistry).PAR)[1, 1, 1:2]

    return all(results_PAR .≈ reverse(expected_PAR))
end

archs = (CPU(), )

@testset "Light attenuaiton model" begin 
    for model in (NonhydrostaticModel, HydrostaticFreeSurfaceModel),
        arch in archs,
        grid in (RectilinearGrid(arch; size = (2, 2, 2), extent = (2, 2, 2)), 
                 LatitudeLongitudeGrid(arch; size = (5, 5, 2), longitude = (-180, 180), latitude = (-85, 85), z = (-2, 0))),
        bgc in (LOBSTER, NutrientPhytoplanktonZooplanktonDetritus)

        if !((model == NonhydrostaticModel) && ((grid isa LatitudeLongitudeGrid) | (grid isa OrthogonalSphericalShellGrid)))
            @info "Testing $bgc in $model on $grid..."
            @test test_two_band(grid, bgc, model)
        end
    end
end
             

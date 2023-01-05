"""
Nutrient-Phytoplankton-Zooplankton-Detritus model of [Kuhn2015](@cite)
http://dx.doi.org/10.1016/j.pocean.2015.07.004
"""

module NPZDModel

export NutrientPhytoplanktonZooplanktonDetritus

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Units
using Oceananigans.Advection: CenteredSecondOrder

using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRatiation, update_PAR!, required_PAR_fields
using OceanBioME: setup_velocity_fields
using OceanBioME.BoxModels: BoxModel

import OceanBioME.BoxModels: update_boxmodel_state!

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    update_biogeochemical_state!,
    biogeochemical_drift_velocity,
    biogeochemical_advection_scheme

struct NutrientPhytoplanktonZooplanktonDetritus{FT, LA, SPAR, W, A} <: AbstractContinuousFormBiogeochemistry
    # phytoplankton
    initial_photosynthetic_slope :: FT # α, 1/(W/m²)/s
    base_maximum_growth :: FT # μ₀, 1/s
    nutrient_half_saturation :: FT # kₙ, mmol N/m³
    base_respiration_rate :: FT # lᵖⁿ, 1/s
    phyto_base_mortality_rate :: FT # lᵖᵈ, 1/s

    # zooplankton
    maximum_grazing_rate :: FT # gₘₐₓ, 1/s
    grazing_half_saturation :: FT # kₚ, mmol N/m³
    assimulation_efficiency :: FT # β
    base_excretion_rate :: FT # lᶻⁿ, 1/s
    zoo_base_mortality_rate :: FT # lᶻᵈ, 1/s

    # detritus
    remineralization_rate :: FT # rᵈⁿ, 1/s

    # light attenuation
    light_attenuation_model :: LA
    surface_phytosynthetically_active_radiation :: SPAR

    # sinking
    sinking_velocities :: W
    advection_schemes :: A

    function NutrientPhytoplanktonZooplanktonDetritus(;grid,
                                                       initial_photosynthetic_slope::FT = 0.1953 / day, # 1/(W/m²)/s
                                                       base_maximum_growth::FT = 0.6989 / day, # 1/s
                                                       nutrient_half_saturation::FT = 2.3868, # mmol N/m³
                                                       base_respiration_rate::FT = 0.066 / day, # 1/s/(mmol N / m³)
                                                       phyto_base_mortality_rate::FT = 0.0101 / day, # 1/s/(mmol N / m³)
                                                       maximum_grazing_rate::FT = 2.1522 / day, # 1/s
                                                       grazing_half_saturation::FT = 0.5573, # mmol N/m³
                                                       assimulation_efficiency::FT = 0.9116, 
                                                       base_excretion_rate::FT = 0.0102 / day, # 1/s/(mmol N / m³)
                                                       zoo_base_mortality_rate::FT = 0.3395 / day, # 1/s/(mmol N / m³)²
                                                       remineralization_rate::FT = 0.1213 / day, # 1/s

                                                       light_attenuation_model::LA = TwoBandPhotosyntheticallyActiveRatiation(),
                                                       surface_phytosynthetically_active_radiation::SPAR = (x, y, t) -> 100*max(0.0, cos(t*π/(12hours))),
                
                                                       sinking_speeds = (P = 0.2551/day, D = 2.7489/day),
                                                       open_bottom::Bool = true,
                                                       advection_schemes::A = NamedTuple{keys(sinking_speeds)}(repeat([CenteredSecondOrder()], 
                                                                                              length(sinking_speeds)))) where {FT, LA, SPAR, A}
        sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)
        W = typeof(sinking_velocities)
        return new{FT, LA, SPAR, W, A}(initial_photosynthetic_slope,
                                       base_maximum_growth,
                                       nutrient_half_saturation,
                                       base_respiration_rate,
                                       phyto_base_mortality_rate,
                                       maximum_grazing_rate,
                                       grazing_half_saturation,
                                       assimulation_efficiency,
                                       base_excretion_rate,
                                       zoo_base_mortality_rate,
                                       remineralization_rate,
                                       light_attenuation_model,
                                       surface_phytosynthetically_active_radiation,
                                       sinking_velocities,
                                       advection_schemes)
    end
end

required_biogeochemical_tracers(::NutrientPhytoplanktonZooplanktonDetritus) = (:N, :P, :Z, :D, :T)
required_biogeochemical_auxiliary_fields(::NutrientPhytoplanktonZooplanktonDetritus) = (:PAR, )

@inline nutrient_limitation(N, kₙ) = N / (kₙ + N)

@inline Q₁₀(T) = 1.88 ^ (T / 10) # T in °C

@inline light_limitation(PAR, T, α, μ₀) = α * PAR / sqrt((μ₀ * Q₁₀(T)) ^ 2 + α ^ 2 * PAR ^ 2)

@inline function (bgc::NutrientPhytoplanktonZooplanktonDetritus)(::Val{:N}, x, y, z, t, N, P, Z, D, T, PAR)
    μ₀ = bgc.base_maximum_growth
    kₙ = bgc.nutrient_half_saturation
    α = bgc.initial_photosynthetic_slope
    lᵖⁿ = bgc.base_respiration_rate
    lᶻⁿ = bgc.base_excretion_rate
    rᵈⁿ = bgc.remineralization_rate

    phytoplankton_consumption = μ₀ * Q₁₀(T) * nutrient_limitation(N, kₙ) * light_limitation(PAR, T, α, μ₀ * Q₁₀(T)) * P
    phytoplankton_metabolic_loss = lᵖⁿ * Q₁₀(T) * P
    zooplankton_metabolic_loss = lᶻⁿ * Q₁₀(T) * Z
    remineralization = rᵈⁿ * D

    return phytoplankton_metabolic_loss + zooplankton_metabolic_loss + remineralization - phytoplankton_consumption
end

@inline function (bgc::NutrientPhytoplanktonZooplanktonDetritus)(::Val{:P}, x, y, z, t, N, P, Z, D, T, PAR)
    μ₀ = bgc.base_maximum_growth
    kₙ = bgc.nutrient_half_saturation
    α = bgc.initial_photosynthetic_slope
    gₘₐₓ = bgc.maximum_grazing_rate
    kₚ = bgc.grazing_half_saturation
    lᵖⁿ = bgc.base_respiration_rate
    lᵖᵈ = bgc.phyto_base_mortality_rate

    growth = μ₀ * Q₁₀(T) * nutrient_limitation(N, kₙ) * light_limitation(PAR, T, α, μ₀ * Q₁₀(T)) * P
    grazing = gₘₐₓ * nutrient_limitation(P ^ 2, kₚ ^ 2) * Z
    metabolic_loss = lᵖⁿ * Q₁₀(T) * P
    mortality_loss = lᵖᵈ * Q₁₀(T) * P

    return growth - grazing - metabolic_loss - mortality_loss
end

@inline function (bgc::NutrientPhytoplanktonZooplanktonDetritus)(::Val{:Z}, x, y, z, t, N, P, Z, D, T, PAR)
    gₘₐₓ = bgc.maximum_grazing_rate
    kₚ = bgc.grazing_half_saturation
    lᶻⁿ = bgc.base_excretion_rate
    lᶻᵈ = bgc.zoo_base_mortality_rate
    β = bgc.assimulation_efficiency

    grazing = β * gₘₐₓ * nutrient_limitation(P ^ 2, kₚ ^ 2) * Z
    metabolic_loss = lᶻⁿ * Q₁₀(T) * Z
    mortality_loss = lᶻᵈ * Q₁₀(T) * Z ^ 2

    return grazing - metabolic_loss - mortality_loss 
end

@inline function (bgc::NutrientPhytoplanktonZooplanktonDetritus)(::Val{:D}, x, y, z, t, N, P, Z, D, T, PAR)
    lᵖᵈ = bgc.phyto_base_mortality_rate
    gₘₐₓ = bgc.maximum_grazing_rate
    kₚ = bgc.grazing_half_saturation
    β = bgc.assimulation_efficiency
    lᶻᵈ = bgc.zoo_base_mortality_rate
    rᵈⁿ = bgc.remineralization_rate

    phytoplankton_mortality_loss = lᵖᵈ * Q₁₀(T) * P
    zooplankton_assimilation_loss = (1 - β) * gₘₐₓ * nutrient_limitation(P ^ 2, kₚ ^ 2) * Z
    zooplankton_mortality_loss = lᶻᵈ * Q₁₀(T) * Z ^ 2
    remineralization = rᵈⁿ * D

    return phytoplankton_mortality_loss + zooplankton_assimilation_loss + zooplankton_mortality_loss - remineralization
end

@inline function biogeochemical_drift_velocity(bgc::NutrientPhytoplanktonZooplanktonDetritus, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return bgc.sinking_velocities[tracer_name]
    else
        return nothing
    end
end

@inline function biogeochemical_advection_scheme(bgc::NutrientPhytoplanktonZooplanktonDetritus, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return bgc.advection_schemes[tracer_name]
    else
        return nothing
    end
end

function update_biogeochemical_state!(bgc::NutrientPhytoplanktonZooplanktonDetritus, model)
    update_PAR!(model, bgc.light_attenuation_model, bgc.surface_phytosynthetically_active_radiation)
end

function update_boxmodel_state!(model::BoxModel{<:NutrientPhytoplanktonZooplanktonDetritus, <:Any, <:Any, <:Any, <:Any, <:Any})
    getproperty(model.values, :PAR) .= model.forcing.PAR(model.clock.time)
    getproperty(model.values, :T) .= model.forcing.T(model.clock.time)
end

end # module
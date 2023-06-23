"""
Nutrient-Phytoplankton-Zooplankton-Detritus model of [Kuhn2015](@cite)

Tracers
========
* Nutrients: N (mmol N/m³)
* Phytoplankton: P (mmol N/m³)
* Zooplankton: Z (mmol N/m³)
* Detritus: D (mmol N/m³)

Required submodels
===========
* Photosynthetically available radiation: PAR (W/m²)

"""
module NPZDModel

export NutrientPhytoplanktonZooplanktonDetritus, NPZD

using OceanBioME: ContinuousFormBiogeochemistry

using Oceananigans.Units

using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRatiation, update_PAR!, required_PAR_fields
using OceanBioME: setup_velocity_fields, show_sinking_velocities
using OceanBioME.BoxModels: BoxModel

import OceanBioME.BoxModels: update_boxmodel_state!
import Base: show, summary

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     update_biogeochemical_state!,
                                     biogeochemical_drift_velocity,
				                     biogeochemical_auxiliary_fields

import OceanBioME: maximum_sinking_velocity

import Adapt: adapt_structure, adapt

"""
    NutrientPhytoplanktonZooplanktonDetritus(; grid,
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

        surface_phytosynthetically_active_radiation = (x, y, t) -> 100*max(0.0, cos(t*π/(12hours))),
        light_attenuation_model::LA = TwoBandPhotosyntheticallyActiveRatiation(; grid,
                                    surface_PAR = surface_phytosynthetically_active_radiation),
        sediment_model::S = nothing,

        sinking_speeds = (P = 0.2551/day, D = 2.7489/day),
        open_bottom::Bool = true,
                                           
        particles::P = nothing)

Construct an instance of the NPZD ([NPZD](@ref NPZD)) biogeochemical model.

Keywork Arguments
===================

- `grid`: (required) the geometry to build the model on, required to calculate sinking
- `initial_photosynthetic_slope`, ..., `remineralization_rate`: NPZD parameter values
- `surface_phytosynthetically_active_radiation`: funciton (or array in the future) for the photosynthetically available radiaiton at the surface, should be shape `f(x, y, t)`
- `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `AbstractSediment`
- `sinking_speed`: named tuple of constant sinking, of fields (i.e. `ZFaceField(...)`) for any tracers which sink (convention is that a sinking speed is positive, but a field will need to follow the usual down being negative)
- `open_bottom`: should the sinking velocity be smoothly brought to zero at the bottom to prevent the tracers leaving the domain
- `particles`: slot for `BiogeochemicalParticles`
"""

struct NutrientPhytoplanktonZooplanktonDetritus{FT, LA, S, W, P} <: ContinuousFormBiogeochemistry{LA, S, P}
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

    # sediment
    sediment_model :: S

    # sinking
    sinking_velocities :: W

    particles :: P

    function NutrientPhytoplanktonZooplanktonDetritus(initial_photosynthetic_slope::FT,
                                                      base_maximum_growth::FT,
                                                      nutrient_half_saturation::FT,
                                                      base_respiration_rate::FT,
                                                      phyto_base_mortality_rate::FT,
    
                                                      maximum_grazing_rate::FT,
                                                      grazing_half_saturation::FT,
                                                      assimulation_efficiency::FT,
                                                      base_excretion_rate::FT,
                                                      zoo_base_mortality_rate::FT,
    
                                                      remineralization_rate::FT,
    
                                                      light_attenuation_model::LA,
    
                                                      sediment_model::S,
    
                                                      sinking_velocities::W,
    
                                                      particles::P) where {FT, LA, S, W, P}
        return new{FT, LA, S, W, P}(initial_photosynthetic_slope,
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

                                    sediment_model,

                                    sinking_velocities,

                                    particles)
    end
end

function NutrientPhytoplanktonZooplanktonDetritus(; grid,
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

                                                    surface_phytosynthetically_active_radiation = (x, y, t) -> 100*max(0.0, cos(t*π/(12hours))),
                                                    light_attenuation_model::LA = TwoBandPhotosyntheticallyActiveRatiation(; grid,
                                                                                    surface_PAR = surface_phytosynthetically_active_radiation),
                                                    sediment_model::S = nothing,
                
                                                    sinking_speeds = (P = 0.2551/day, D = 2.7489/day),
                                                    open_bottom::Bool = true,
                                                                                           
                                                    particles::P = nothing) where {FT, LA, S, A, P}

    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)

    W = typeof(sinking_velocities)

    return NutrientPhytoplanktonZooplanktonDetritus(initial_photosynthetic_slope,
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
                                                    sediment_model,
                                                    sinking_velocities,
                                                    particles)
end

const NPZD = NutrientPhytoplanktonZooplanktonDetritus

required_biogeochemical_tracers(::NPZD) = (:N, :P, :Z, :D, :T)
required_biogeochemical_auxiliary_fields(::NPZD) = (:PAR, )

@inline nutrient_limitation(N, kₙ) = N / (kₙ + N)

@inline Q₁₀(T) = 1.88 ^ (T / 10) # T in °C

@inline light_limitation(PAR, T, α, μ₀) = α * PAR / sqrt((μ₀ * Q₁₀(T)) ^ 2 + α ^ 2 * PAR ^ 2)

@inline function (bgc::NPZD)(::Val{:N}, x, y, z, t, N, P, Z, D, T, PAR)
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

@inline function (bgc::NPZD)(::Val{:P}, x, y, z, t, N, P, Z, D, T, PAR)
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

@inline function (bgc::NPZD)(::Val{:Z}, x, y, z, t, N, P, Z, D, T, PAR)
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

@inline function (bgc::NPZD)(::Val{:D}, x, y, z, t, N, P, Z, D, T, PAR)
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

@inline function biogeochemical_drift_velocity(bgc::NPZD, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return bgc.sinking_velocities[tracer_name]
    else
        return nothing
    end
end

function update_biogeochemical_state!(bgc::NPZD, model)
    update_PAR!(model, bgc.light_attenuation_model)
end

function update_boxmodel_state!(model::BoxModel{<:NPZD, <:Any, <:Any, <:Any, <:Any, <:Any})
    getproperty(model.values, :PAR) .= model.forcing.PAR(model.clock.time)
    getproperty(model.values, :T) .= model.forcing.T(model.clock.time)
end

summary(::NPZD{FT, LA, W, A, P}) where {FT, LA, W, A, P} = string("Nutrient Phytoplankton Zooplankton Detritus model ($FT)")
show(io::IO, model::NPZD{FT, LA, W, A, P}) where {FT, LA, W, A, P} =
       print(io, summary(model), " \n",
                " Light Attenuation Model: ", "\n",
                "    └── ", summary(model.light_attenuation_model), "\n",
                " Sinking Velocities:", "\n", show_sinking_velocities(model.sinking_velocities))

@inline maximum_sinking_velocity(bgc::NPZD) = maximum(abs, bgc.sinking_velocities.D.w)
@inline biogeochemical_auxiliary_fields(bgc::NPZD) = biogeochemical_auxiliary_fields(bgc.light_attenuation_model)

adapt_structure(to, npzd::NPZD) = 
    NutrientPhytoplanktonZooplanktonDetritus(npzd.initial_photosynthetic_slope,
                                             npzd.base_maximum_growth,
                                             npzd.nutrient_half_saturation,
                                             npzd.base_respiration_rate,
                                             npzd.phyto_base_mortality_rate,

                                             npzd.maximum_grazing_rate,
                                             npzd.grazing_half_saturation,
                                             npzd.assimulation_efficiency,
                                             npzd.base_excretion_rate,
                                             npzd.zoo_base_mortality_rate,

                                             npzd.remineralization_rate,

                                             adapt(to, npzd.light_attenuation_model),
                                             adapt(to, npzd.sediment_model),

                                             adapt(to, npzd.sinking_velocities),

                                             adapt(to, npzd.particles))
end # module

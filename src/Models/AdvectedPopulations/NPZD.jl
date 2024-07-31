"""
Nutrient-Phytoplankton-Zooplankton-Detritus model of [Kuhn2015](@citet).

Tracers
=======
* Nutrients: N (mmol N/m³)
* Phytoplankton: P (mmol N/m³)
* Zooplankton: Z (mmol N/m³)
* Detritus: D (mmol N/m³)

Required submodels
==================
* Photosynthetically available radiation: PAR (W/m²)
"""
module NPZDModel

export NutrientPhytoplanktonZooplanktonDetritus, NPZD

using OceanBioME: Biogeochemistry, ScaleNegativeTracers
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry

using Oceananigans.Units
using Oceananigans.Fields: ZeroField

using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRadiation, default_surface_PAR
using OceanBioME: setup_velocity_fields, show_sinking_velocities
using OceanBioME.BoxModels: BoxModel
using OceanBioME.Boundaries.Sediments: sinking_flux

import Base: show, summary

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers,
                                     required_biogeochemical_auxiliary_fields,
                                     biogeochemical_drift_velocity

import OceanBioME: redfield, conserved_tracers, maximum_sinking_velocity, chlorophyll
import OceanBioME.Boundaries.Sediments: nitrogen_flux, carbon_flux, remineralisation_receiver, sinking_tracers

import Adapt: adapt_structure, adapt

struct NutrientPhytoplanktonZooplanktonDetritus{FT, W} <: AbstractContinuousFormBiogeochemistry
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

    # sinking
    sinking_velocities :: W

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
    
                                                      sinking_velocities::W) where {FT, W}
        return new{FT, W}(initial_photosynthetic_slope,
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
                          
                          sinking_velocities)
    end
end

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

                                               surface_photosynthetically_active_radiation = default_surface_PAR,
                                               light_attenuation_model::LA =
                                                   TwoBandPhotosyntheticallyActiveRadiation(; grid,
                                                                                              surface_PAR = surface_photosynthetically_active_radiation),
                                               sediment_model::S = nothing,
                
                                               sinking_speeds = (P = 0.2551/day, D = 2.7489/day),
                                               open_bottom::Bool = true,

                                               scale_negatives = false,
                                                                                      
                                               particles::P = nothing,
                                               modifiers::M = nothing)

Construct a Nutrient-Phytoplankton-Zooplankton-Detritus ([NPZD](@ref NPZD)) biogeochemical model.

Keyword Arguments
=================

- `grid`: (required) the geometry to build the model on, required to calculate sinking
- `initial_photosynthetic_slope`, ..., `remineralization_rate`: NPZD parameter values
- `surface_photosynthetically_active_radiation`: function (or array in the future) for the photosynthetically available radiation at the surface, should be shape `f(x, y, t)`
- `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
- `sediment_model`: slot for `AbstractSediment`
- `sinking_speed`: named tuple of constant sinking, of fields (i.e. `ZFaceField(...)`) for any tracers which sink (convention is that a sinking speed is positive, but a field will need to follow the usual down being negative)
- `open_bottom`: should the sinking velocity be smoothly brought to zero at the bottom to prevent the tracers leaving the domain
- `scale_negatives`: scale negative tracers?
- `particles`: slot for `BiogeochemicalParticles`
- `modifiers`: slot for components which modify the biogeochemistry when the tendencies have been calculated or when the state is updated

Example
=======

```jldoctest
julia> using OceanBioME

julia> using Oceananigans

julia> grid = RectilinearGrid(size=(20, 30), extent=(200, 200), topology=(Bounded, Flat, Bounded));

julia> model = NutrientPhytoplanktonZooplanktonDetritus(; grid)
NutrientPhytoplanktonZooplanktonDetritus{Float64} model, with (:P, :D) sinking 
 Light attenuation: Two-band light attenuation model (Float64)
 Sediment: Nothing
 Particles: Nothing
 Modifiers: Nothing
```
"""
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

                                                    surface_photosynthetically_active_radiation = default_surface_PAR,
                                                    light_attenuation_model::LA =
                                                        TwoBandPhotosyntheticallyActiveRadiation(; grid,
                                                                                                   surface_PAR = surface_photosynthetically_active_radiation),
                                                    sediment_model::S = nothing,
                
                                                    sinking_speeds = (P = 0.2551/day, D = 2.7489/day),
                                                    open_bottom::Bool = true,

                                                    scale_negatives = false,
                                                    invalid_fill_value = NaN,
                                                                                           
                                                    particles::P = nothing,
                                                    modifiers::M = nothing) where {FT, LA, S, P, M}

    sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)

    underlying_biogeochemistry = 
        NutrientPhytoplanktonZooplanktonDetritus(initial_photosynthetic_slope,
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
                                                 sinking_velocities)

    if scale_negatives
        scaler = ScaleNegativeTracers(underlying_biogeochemistry, grid; invalid_fill_value)
        modifiers = isnothing(modifiers) ? scaler : (modifiers..., scaler)
    end

    return Biogeochemistry(underlying_biogeochemistry;
                           light_attenuation = light_attenuation_model, 
                           sediment = sediment_model, 
                           particles,
                           modifiers)
end

const NPZD = NutrientPhytoplanktonZooplanktonDetritus

required_biogeochemical_tracers(::NPZD) = (:N, :P, :Z, :D, :T)
required_biogeochemical_auxiliary_fields(::NPZD) = (:PAR, )

@inline nutrient_limitation(N, kₙ) = N / (kₙ + N)

@inline Q₁₀(T) = 1.88 ^ (T / 10) # T in °C

@inline light_limitation(PAR, α, μ₀) = α * PAR / sqrt(μ₀ ^ 2 + α ^ 2 * PAR ^ 2)

@inline function (bgc::NPZD)(::Val{:N}, x, y, z, t, N, P, Z, D, T, PAR)
    μ₀ = bgc.base_maximum_growth
    kₙ = bgc.nutrient_half_saturation
    α = bgc.initial_photosynthetic_slope
    lᵖⁿ = bgc.base_respiration_rate
    lᶻⁿ = bgc.base_excretion_rate
    rᵈⁿ = bgc.remineralization_rate

    phytoplankton_consumption = μ₀ * Q₁₀(T) * nutrient_limitation(N, kₙ) * light_limitation(PAR, α, μ₀ * Q₁₀(T)) * P
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

    growth = μ₀ * Q₁₀(T) * nutrient_limitation(N, kₙ) * light_limitation(PAR, α, μ₀ * Q₁₀(T)) * P
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
        return (u = ZeroField(), v = ZeroField(), w = bgc.sinking_velocities[tracer_name])
    else
        return (u = ZeroField(), v = ZeroField(), w = ZeroField())
    end
end

summary(::NPZD{FT, NamedTuple{K, V}}) where {FT, K, V} = string("NutrientPhytoplanktonZooplanktonDetritus{$FT} model, with $K sinking")
show(io::IO, model::NPZD{FT}) where {FT} = print(io, string("NutrientPhytoplanktonZooplanktonDetritus{$FT} model \n",
                                                            "└── Sinking Velocities:", "\n", show_sinking_velocities(model.sinking_velocities)))

@inline maximum_sinking_velocity(bgc::NPZD) = maximum(abs, bgc.sinking_velocities.D.w)

adapt_structure(to, npzd::NPZD) = 
    NutrientPhytoplanktonZooplanktonDetritus(adapt(to, npzd.initial_photosynthetic_slope),
                                             adapt(to, npzd.base_maximum_growth),
                                             adapt(to, npzd.nutrient_half_saturation),
                                             adapt(to, npzd.base_respiration_rate),
                                             adapt(to, npzd.phyto_base_mortality_rate),

                                             adapt(to, npzd.maximum_grazing_rate),
                                             adapt(to, npzd.grazing_half_saturation),
                                             adapt(to, npzd.assimulation_efficiency),
                                             adapt(to, npzd.base_excretion_rate),
                                             adapt(to, npzd.zoo_base_mortality_rate),

                                             adapt(to, npzd.remineralization_rate),

                                             adapt(to, npzd.sinking_velocities))

@inline redfield(i, j, k, val_tracer_name, bgc::NPZD, tracers) = redfield(val_tracer_name, bgc)
@inline redfield(::Union{Val{:N}}, bgc::NPZD) = 0
@inline redfield(::Union{Val{:P}, Val{:Z}, Val{:D}}, bgc::NPZD) = 6.56

@inline nitrogen_flux(i, j, k, grid, advection, bgc::NPZD, tracers) = sinking_flux(i, j, k, grid, advection, Val(:D), bgc, tracers) +
                                                                      sinking_flux(i, j, k, grid, advection, Val(:P), bgc, tracers)

@inline carbon_flux(i, j, k, grid, advection, bgc::NPZD, tracers) = nitrogen_flux(i, j, k, grid, advection, bgc, tracers) * redfield(Val(:P), bgc)

@inline remineralisation_receiver(::NPZD) = :N

@inline conserved_tracers(::NPZD) = (:N, :P, :Z, :D)
@inline sinking_tracers(bgc::NPZD) = keys(bgc.sinking_velocities)

@inline chlorophyll(bcg::NPZD, model) = 1.31 * model.tracers.P

end # module

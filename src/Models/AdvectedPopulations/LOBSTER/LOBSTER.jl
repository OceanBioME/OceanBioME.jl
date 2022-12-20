"""
The Lodyc Ocean Biogeochemical Simulation Tools for Ecosystem and Resources (LOBSTER) model

Tracers
========
* Nitrates: NO₃ (mmol N/m³)
* Ammonia: NH₄ (mmol N/m³)
* Phytoplankton: P (mmol N/m³)
* Zooplankton: Z (mmol N/m³)
* Small (slow sinking) particulate organic matter: sPOM (mmol N/m³)
* Large (fast sinking) particulate organic matter: bPOM (mmol N/m³)
* Disolved organic matter: DOM (mmol N/m³)

Optional tracers
===========
Carbonate chemistry
* Disolved inorganic carbon: DIC (mmol C/m³)
* Alkalinity: Alk (meq/m³)

Oxygen chemistry
* Oxygen: O₂ (mmol O₂/m³)

Variable redfield
* Small (slow sinking) particulate organic matter carbon content: sPOC (mmol C/m³)
* Large (fast sinking) particulate organic matter carbon content: bPOC (mmol C/m³)
* Disolved organic matter carbon content: DOC (mmol C/m³)
* When this option is enabled then the usual sPOM and bPOM change to sPON and bPON as they explicitly represent the nitrogen contained in the particulate matter

Required forcing
===========
* Photosynthetically available radiation: PAR (W/m²)

For optional tracers:
* Temperature: T (ᵒC)
* Salinity: S (‰)
"""
module LOBSTERModel

export LOBSTER

using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
using Oceananigans.Units
using Oceananigans.Advection: CenteredSecondOrder
using Oceananigans.Fields: Field, TracerFields, CenterField

using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRatiation, update_PAR!, required_PAR_fields
using OceanBioME: setup_velocity_fields
using OceanBioME.BoxModels: BoxModel

import OceanBioME.BoxModels: update_boxmodel_state!

import Oceananigans.Biogeochemistry:
       required_biogeochemical_tracers,
       required_biogeochemical_auxiliary_fields,
       biogeochemical_drift_velocity,
       biogeochemical_advection_scheme,
       update_biogeochemical_state!

"""
    LOBSTER(;grid,
             phytoplankton_preference::FT = 0.5,
             maximum_grazing_rate::FT = 9.26e-6, # 1/s
             grazing_half_saturation::FT = 1.0, # mmol N/m³
             light_half_saturation::FT = 33.0, # W/m² (?)
             nitrate_ammonia_inhibition::FT = 3.0,
             nitrate_half_saturation::FT = 0.7, # mmol N/m³
             ammonia_half_saturation::FT = 0.001, # mmol N/m³
             maximum_phytoplankton_growthrate::FT = 1.21e-5, # 1/s
             zooplankton_assimilation_fraction::FT = 0.7,
             zooplankton_mortality::FT = 2.31e-6, # 1/s/mmol N/m³
             zooplankton_excretion_rate::FT = 5.8e-7, # 1/s
             phytoplankton_mortality::FT = 5.8e-7, # 1/s
             small_detritus_remineralisation_rate::FT = 5.88e-7, # 1/s
             large_detritus_remineralisation_rate::FT = 5.88e-7, # 1/s
             phytoplankton_exudation_fraction::FT = 0.05,
             nitrifcaiton_rate::FT = 5.8e-7, # 1/s
             ammonia_fraction_of_exudate::FT = 0.75, 
             ammonia_fraction_of_excriment::FT = 0.5,
             ammonia_fraction_of_detritus::FT = 0.0,
             phytoplankton_redfield::FT = 6.56, # mol C/mol N
             organic_redfield::FT = 6.56, # mol C/mol N
             phytoplankton_chlorophyll_ratio::FT = 1.31, # g Chl/mol N
             organic_carbon_calcate_ratio::FT = 0.1, # mol CaCO₃/mol C
             respiraiton_oxygen_nitrogen_ratio::FT = 10.75, # mol O/molN
             nitrifcation_oxygen_nitrogen_ratio::FT = 2.0, # mol O/molN
             slow_sinking_mortality_fraction::FT = 0.5, 
             fast_sinking_mortality_fraction::FT = 0.5,
             disolved_organic_breakdown_rate::FT = 3.86e-7, # 1/s
             zooplankton_calcite_dissolution::FT = 0.3,

             light_attenuation_model = TwoBandPhotosyntheticallyActiveRatiation(),
             surface_phytosynthetically_active_radiation::SPAR = (x, y, t) -> 100*max(0.0, cos(t*π/(12hours))),

             carbonates::Bool = false,
             oxygen::Bool = false,
             variable_redfield = false,

             sinking_speed = (sPOM = 3.47e-5, bPOM = 200/day),
             open_bottom::Bool = true,
             advection_schemes::A = NamedTuple{keys(sinking_velocities)}(repeat([CenteredSecondOrder()], 
                                                                         length(sinking_velocities))))

Construct an instance of the LOBSTER ([LOBSTER](@ref LOBSTER)) biogeochemical model.

Keywork Arguments
===================

    - `grid`: (required) the geometry to build the model on, required to calculate sinking
    - `phytoplankton_preference`, ..., `disolved_organic_breakdown_rate`: LOBSTER parameter values
    - `light_attenuation_model`: light attenuation model which integrated the attenuation of available light
    - `surface_phytosynthetically_active_radiation`: funciton (or array in the future) for the photosynthetically available radiaiton at the surface, should be shape `f(x, y, t)`
    - `carbonates`, `oxygen`, and `variable_redfield`: include models for carbonate chemistry and/or oxygen chemistry and/or variable redfield ratio disolved and particulate organic matter
    - `sinking_speed`: named tuple of constant sinking, of fields (i.e. `ZFaceField(...)`) for any tracers which sink (convention is that a sinking speed is positive, but a field will need to follow the usual down being negative)
    - `open_bottom`: should the sinking velocity be smoothly brought to zero at the bottom to prevent the tracers leaving the domain
    - `advection_schemes`: named tuple of advection scheme to use for sinking
"""
struct LOBSTER{FT, LA, SPAR, B, W, A} <: AbstractContinuousFormBiogeochemistry
    phytoplankton_preference :: FT
    maximum_grazing_rate :: FT
    grazing_half_saturation :: FT
    light_half_saturation :: FT
    nitrate_ammonia_inhibition :: FT
    nitrate_half_saturation :: FT
    ammonia_half_saturation :: FT
    maximum_phytoplankton_growthrate :: FT
    zooplankton_assimilation_fraction :: FT
    zooplankton_mortality :: FT
    zooplankton_excretion_rate :: FT
    phytoplankton_mortality :: FT
    small_detritus_remineralisation_rate :: FT
    large_detritus_remineralisation_rate :: FT
    phytoplankton_exudation_fraction :: FT
    nitrifcaiton_rate :: FT
    ammonia_fraction_of_exudate :: FT
    ammonia_fraction_of_excriment :: FT
    ammonia_fraction_of_detritus :: FT
    phytoplankton_redfield :: FT
    organic_redfield :: FT
    phytoplankton_chlorophyll_ratio :: FT
    organic_carbon_calcate_ratio :: FT
    respiraiton_oxygen_nitrogen_ratio :: FT
    nitrifcation_oxygen_nitrogen_ratio :: FT
    slow_sinking_mortality_fraction :: FT
    fast_sinking_mortality_fraction :: FT
    disolved_organic_breakdown_rate :: FT
    zooplankton_calcite_dissolution :: FT

    light_attenuation_model :: LA
    surface_phytosynthetically_active_radiation :: SPAR

    optionals :: B

    sinking_velocities :: W
    advection_schemes :: A

    function LOBSTER(;grid,
                      phytoplankton_preference::FT = 0.5,
                      maximum_grazing_rate::FT = 9.26e-6, # 1/s
                      grazing_half_saturation::FT = 1.0, # mmol N/m³
                      light_half_saturation::FT = 33.0, # W/m² (?)
                      nitrate_ammonia_inhibition::FT = 3.0,
                      nitrate_half_saturation::FT = 0.7, # mmol N/m³
                      ammonia_half_saturation::FT = 0.001, # mmol N/m³
                      maximum_phytoplankton_growthrate::FT = 1.21e-5, # 1/s
                      zooplankton_assimilation_fraction::FT = 0.7,
                      zooplankton_mortality::FT = 2.31e-6, # 1/s/mmol N/m³
                      zooplankton_excretion_rate::FT = 5.8e-7, # 1/s
                      phytoplankton_mortality::FT = 5.8e-7, # 1/s
                      small_detritus_remineralisation_rate::FT = 5.88e-7, # 1/s
                      large_detritus_remineralisation_rate::FT = 5.88e-7, # 1/s
                      phytoplankton_exudation_fraction::FT = 0.05,
                      nitrifcaiton_rate::FT = 5.8e-7, # 1/s
                      ammonia_fraction_of_exudate::FT = 0.75, 
                      ammonia_fraction_of_excriment::FT = 0.5,
                      ammonia_fraction_of_detritus::FT = 0.0,
                      phytoplankton_redfield::FT = 6.56, # mol C/mol N
                      organic_redfield::FT = 6.56, # mol C/mol N
                      phytoplankton_chlorophyll_ratio::FT = 1.31, # g Chl/mol N
                      organic_carbon_calcate_ratio::FT = 0.1, # mol CaCO₃/mol C
                      respiraiton_oxygen_nitrogen_ratio::FT = 10.75, # mol O/molN
                      nitrifcation_oxygen_nitrogen_ratio::FT = 2.0, # mol O/molN
                      slow_sinking_mortality_fraction::FT = 0.5, 
                      fast_sinking_mortality_fraction::FT = 0.5,
                      disolved_organic_breakdown_rate::FT = 3.86e-7, # 1/s
                      zooplankton_calcite_dissolution::FT = 0.3,

                      light_attenuation_model::LA = TwoBandPhotosyntheticallyActiveRatiation(),
                      surface_phytosynthetically_active_radiation::SPAR = (x, y, t) -> 100*max(0.0, cos(t*π/(12hours))),

                      carbonates::Bool = false,
                      oxygen::Bool = false,
                      variable_redfield::Bool = false,
                
                      sinking_speeds = (sPOM = 3.47e-5, bPOM = 200/day),
                      open_bottom::Bool = true,
                      advection_schemes::A = NamedTuple{keys(sinking_speeds)}(repeat([CenteredSecondOrder(grid)], 
                                                                                    length(sinking_speeds)))) where {FT, LA, SPAR, A}

        

        sinking_velocities = setup_velocity_fields(sinking_speeds, grid, open_bottom)
        W = typeof(sinking_velocities)
        optionals = Val((carbonates, oxygen, variable_redfield))
        B = typeof(optionals)

        return new{FT, LA, SPAR, B, W, A}(phytoplankton_preference,
                                          maximum_grazing_rate,
                                          grazing_half_saturation,
                                          light_half_saturation,
                                          nitrate_ammonia_inhibition,
                                          nitrate_half_saturation,
                                          ammonia_half_saturation,
                                          maximum_phytoplankton_growthrate,
                                          zooplankton_assimilation_fraction,
                                          zooplankton_mortality,
                                          zooplankton_excretion_rate,
                                          phytoplankton_mortality,
                                          small_detritus_remineralisation_rate,
                                          large_detritus_remineralisation_rate,
                                          phytoplankton_exudation_fraction,
                                          nitrifcaiton_rate,
                                          ammonia_fraction_of_exudate,
                                          ammonia_fraction_of_excriment,
                                          ammonia_fraction_of_detritus,
                                          phytoplankton_redfield,
                                          organic_redfield,
                                          phytoplankton_chlorophyll_ratio,
                                          organic_carbon_calcate_ratio,
                                          respiraiton_oxygen_nitrogen_ratio,
                                          nitrifcation_oxygen_nitrogen_ratio,
                                          slow_sinking_mortality_fraction,
                                          fast_sinking_mortality_fraction,
                                          disolved_organic_breakdown_rate,
                                          zooplankton_calcite_dissolution,

                                          light_attenuation_model,
                                          surface_phytosynthetically_active_radiation,

                                          optionals,
                                                            
                                          sinking_velocities,
                                          advection_schemes)
    end
end

# wrote this functionally and it took 2.5x longer so even though this is long going to use this way instead
@inline required_biogeochemical_tracers(::LOBSTER{<:Any, <:Any, <:Any, <:Val{(false, false, false)}, <:Any, <:Any}) = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)
@inline required_biogeochemical_tracers(::LOBSTER{<:Any, <:Any, <:Any, <:Val{(true, false, false)}, <:Any, <:Any}) = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM, :DIC, :Alk)
@inline required_biogeochemical_tracers(::LOBSTER{<:Any, <:Any, <:Any, <:Val{(false, true, false)}, <:Any, <:Any}) = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM, :O₂)
@inline required_biogeochemical_tracers(::LOBSTER{<:Any, <:Any, <:Any, <:Val{(false, false, true)}, <:Any, <:Any}) = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON, :sPOC, :bPOC, :DOC)
@inline required_biogeochemical_tracers(::LOBSTER{<:Any, <:Any, <:Any, <:Val{(true, true, false)}, <:Any, <:Any}) = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM, :DIC, :Alk, :O₂)
@inline required_biogeochemical_tracers(::LOBSTER{<:Any, <:Any, <:Any, <:Val{(true, false, true)}, <:Any, <:Any}) = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON, :DIC, :Alk, :sPOC, :bPOC, :DOC)
@inline required_biogeochemical_tracers(::LOBSTER{<:Any, <:Any, <:Any, <:Val{(false, true, true)}, <:Any, <:Any}) = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON, :O₂, :sPOC, :bPOC, :DOC)
@inline required_biogeochemical_tracers(::LOBSTER{<:Any, <:Any, <:Any, <:Val{(true, true, true)}, <:Any, <:Any}) = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON, :DIC, :Alk, :O₂, :sPOC, :bPOC, :DOC)

@inline required_biogeochemical_auxiliary_fields(model::LOBSTER) = required_PAR_fields(model.light_attenuation_model)

const small_detritus = Union{Val{:sPON}, Val{:sPOC}}
const large_detritus = Union{Val{:bPON}, Val{:bPOC}}
const disolved_organic_matter = Union{Val{:DON}, Val{:DOC}}

const sPOM = Union{Val{:sPOM}, Val{:sPON}}
const bPOM = Union{Val{:bPOM}, Val{:bPON}}
const DOM = Union{Val{:DOM}, Val{:DON}}

@inline biogeochemical_drift_velocity(bgc::LOBSTER, ::small_detritus) = biogeochemical_drift_velocity(bgc, Val(:sPOM))
@inline biogeochemical_drift_velocity(bgc::LOBSTER, ::large_detritus) = biogeochemical_drift_velocity(bgc, Val(:bPOM))
@inline biogeochemical_drift_velocity(bgc::LOBSTER, ::disolved_organic_matter) = biogeochemical_drift_velocity(bgc, Val(:DOM))

@inline function biogeochemical_drift_velocity(bgc::LOBSTER, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return bgc.sinking_velocities[tracer_name]
    else
        return nothing
    end
end

@inline biogeochemical_advection_scheme(bgc::LOBSTER, ::small_detritus) = biogeochemical_advection_scheme(bgc, Val(:sPOM))
@inline biogeochemical_advection_scheme(bgc::LOBSTER, ::large_detritus) = biogeochemical_advection_scheme(bgc, Val(:bPOM))
@inline biogeochemical_advection_scheme(bgc::LOBSTER, ::disolved_organic_matter) = biogeochemical_advection_scheme(bgc, Val(:DOM))

@inline function biogeochemical_advection_scheme(bgc::LOBSTER, ::Val{tracer_name}) where tracer_name
    if tracer_name in keys(bgc.sinking_velocities)
        return bgc.advection_schemes[tracer_name]
    else
        return nothing
    end
end

function update_biogeochemical_state!(bgc::LOBSTER, model)
    update_PAR!(model, bgc.light_attenuation_model, bgc.surface_phytosynthetically_active_radiation)
end

function update_boxmodel_state!(model::BoxModel{<:LOBSTER, <:Any, <:Any, <:Any, <:Any, <:Any})
    getproperty(model.values, :PAR) .= model.forcing.PAR(model.clock.time)
end

include("fallbacks.jl")

include("core.jl")
include("carbonate_chemistry.jl")
include("oxygen_chemistry.jl")
include("variable_redfield.jl")

end # module
module OxygenModels

export Oxygen

using ..NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus,
    inorganic_nitrogen_waste

using ..NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    NitrateAmmonia,
    nutrient_uptake,
    nitrification,
    SingleTracerNutrient

using ..NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    primary_production,
    inorganic_carbon_waste

using ..NutrientsPlanktonDetritusModels.DetritusModels: CarbonNitrogenDissolvedParticulate

import Base: summary, show

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields

import ..NutrientsPlanktonDetritusModels: 
    carbon_ratio,
    nitrogen_ratio

"""
    Oxygen([FT = Float64;]
           production_oxygen_carbon_ratio = 131/122,
           nitrification_oxygen_carbon_ratio = 31/122,
           low_oxygen_limit = nothing)

An optional oxygen component for the `oxygen` slot of any [`NutrientsPlanktonDetritus`](@ref) model
(including the [`LOBSTER`](@ref), [`NPZD`](@ref), and [`ImplicitBiology`](@ref)
presets). It adds one tracer, oxygen (`O₂`), which is produced by photosynthesis and consumed by
remineralisation of organic waste and by nitrification.

Oxygen concentration is only one-way coupled to the rest of the biogeochemistry and *does not* affect
any other groups (e.g. low oxygen does *not* reduce zooplankton growth). To capture such an effect a
different `plankton` component would be needed.

Keyword Arguments
=================

- `production_oxygen_carbon_ratio`: moles of O₂ produced per mole of carbon fixed (mol O₂ / mol C)
- `nitrification_oxygen_carbon_ratio`: moles of O₂ consumed per mole of carbon during nitrification
  (mol O₂ / mol C)
- `low_oxygen_limit`: the O₂ concentration (mmol O₂ / m³) at and below which the remineralisation of organic
  waste consumes no O₂; the consumption rises linearly from none there to all of it at twice the limit
  (see [`ramp`](@ref)). The remineralisation itself is not changed, only the O₂ it takes: where O₂ is low
  the organic matter is taken to be oxidised by something else, as by the nitrate of denitrification in a
  model without nitrogen, and the O₂ cannot be driven below the limit by remineralisation. MARBL's is 5 (its
  `parm_o2_min` and `parm_o2_min_delta`). The default, `nothing`, has no limit
"""
struct Oxygen{FT, L}
       production_oxygen_carbon_ratio :: FT
    nitrification_oxygen_carbon_ratio :: FT
                     low_oxygen_limit :: L
end

Oxygen(production_oxygen_carbon_ratio, nitrification_oxygen_carbon_ratio) =
    Oxygen(production_oxygen_carbon_ratio, nitrification_oxygen_carbon_ratio, nothing)

Oxygen(FT::DataType = Float64;
       production_oxygen_carbon_ratio = 131/122,    # mol O₂ / mol C
       nitrification_oxygen_carbon_ratio = 31/122,  # mol O₂ / mol C
       low_oxygen_limit = nothing) =                # mmol O₂ / m³
    Oxygen(convert(FT, production_oxygen_carbon_ratio),
           convert(FT, nitrification_oxygen_carbon_ratio),
           isnothing(low_oxygen_limit) ? nothing : convert(FT, low_oxygen_limit))

"""
    ramp(limit, O₂)

The fraction of the remineralisation's oxygen demand which is taken from the O₂: none at and below `limit`,
all of it at and above twice `limit`, linear between; all of it everywhere when `limit` is `nothing`.
"""
@inline ramp(::Nothing, O₂) = one(O₂)
@inline ramp(limit, O₂) = clamp((O₂ - limit) / limit, zero(O₂), one(O₂))

required_biogeochemical_tracers(::Oxygen) = (:O₂, )
required_biogeochemical_auxiliary_fields(::Oxygen) = tuple()

const OxygenNPD{FT, NUT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:Any, <:Any, <:Oxygen}

# since we're assuming that nitrogen is nitrate (maybe erroneously), then it has to be denitrified before production/after remin
@inline function (bgc::OxygenNPD)(i, j, k, grid, ::Val{:O₂}, clock, fields, auxiliary_fields)
    rP = bgc.oxygen.production_oxygen_carbon_ratio
    rN = bgc.oxygen.nitrification_oxygen_carbon_ratio

    O₂ = @inbounds fields.O₂[i, j, k]
    consumed = ramp(bgc.oxygen.low_oxygen_limit, O₂)

    return (rP + rN) * (
        primary_production(i, j, k, grid,bgc.plankton, bgc, fields, auxiliary_fields)
      - consumed * inorganic_carbon_waste(i, j, k, grid,bgc.plankton, bgc, fields, auxiliary_fields)
      - consumed * inorganic_carbon_waste(i, j, k, grid,bgc.detritus, bgc, fields, auxiliary_fields)
    )
end

@inline function (bgc::OxygenNPD{<:Any, <:Nutrients{<:NitrateAmmonia}})(i, j, k, grid, ::Val{:O₂}, clock, fields, auxiliary_fields)
    rC = carbon_ratio(i, j, k, grid,bgc.plankton, bgc,fields)

    rP = bgc.oxygen.production_oxygen_carbon_ratio
    rN = bgc.oxygen.nitrification_oxygen_carbon_ratio

    # PP from NO₃ - nitrification etc.
    net_nitrate_production = bgc(i, j, k, grid, Val(:NO₃), clock, fields, auxiliary_fields)

    O₂ = @inbounds fields.O₂[i, j, k]
    consumed = ramp(bgc.oxygen.low_oxygen_limit, O₂)

    return (
        rP * primary_production(i, j, k, grid,bgc.plankton, bgc, fields, auxiliary_fields)
      - rN * rC * net_nitrate_production
      - consumed * rP * inorganic_carbon_waste(i, j, k, grid,bgc.plankton, bgc, fields, auxiliary_fields)
      - consumed * rP * inorganic_carbon_waste(i, j, k, grid,bgc.detritus, bgc, fields, auxiliary_fields)
    )
end

const OxygenCNDP{FT} = NutrientsPlanktonDetritus{FT, <:Nutrients{<:SingleTracerNutrient}, <:Any, <:CarbonNitrogenDissolvedParticulate, <:Any, <:Oxygen}

# we will have to think how to generalise this in the future for variable redfield planktons
@inline function (bgc::OxygenCNDP)(i, j, k, grid, ::Val{:O₂}, clock, fields, auxiliary_fields)
    rP = bgc.oxygen.production_oxygen_carbon_ratio
    rN = bgc.oxygen.nitrification_oxygen_carbon_ratio

    O₂ = @inbounds fields.O₂[i, j, k]
    consumed = ramp(bgc.oxygen.low_oxygen_limit, O₂)

    return (
        (rP + rN) * (primary_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
                     - consumed * inorganic_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields))
        - consumed * rP * inorganic_carbon_waste(i, j, k, grid,bgc.detritus, bgc, fields, auxiliary_fields)
        - consumed * rN * inorganic_nitrogen_waste(i, j, k, grid,bgc.detritus, bgc, fields, auxiliary_fields) * carbon_ratio(i, j, k, grid,bgc.plankton, bgc,fields) / nitrogen_ratio(i, j, k, grid,bgc.plankton, bgc,fields)
    )
end
end

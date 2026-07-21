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
           nitrification_oxygen_carbon_ratio = 31/122)

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
"""
struct Oxygen{FT}
       production_oxygen_carbon_ratio :: FT
    nitrification_oxygen_carbon_ratio :: FT
end

Oxygen(FT = Float64;
       production_oxygen_carbon_ratio = 131/122,     # mol O₂ / mol C
       nitrification_oxygen_carbon_ratio = 31/122) = # mol O₂ / mol C
    Oxygen(convert(FT, production_oxygen_carbon_ratio),
           convert(FT, nitrification_oxygen_carbon_ratio))

required_biogeochemical_tracers(::Oxygen) = (:O₂, )
required_biogeochemical_auxiliary_fields(::Oxygen) = tuple()

const OxygenNPD{FT, NUT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:Any, <:Any, <:Oxygen}

# since we're assuming that nitrogen is nitrate (maybe erroneously), then it has to be denitrified before production/after remin
@inline function (bgc::OxygenNPD)(i, j, k, grid, ::Val{:O₂}, clock, fields, auxiliary_fields)
    rP = bgc.oxygen.production_oxygen_carbon_ratio
    rN = bgc.oxygen.nitrification_oxygen_carbon_ratio

    return (rP + rN) * (
        primary_production(i, j, k, grid,bgc.plankton, bgc, fields, auxiliary_fields)
      - inorganic_carbon_waste(i, j, k, grid,bgc.plankton, bgc, fields, auxiliary_fields)
      - inorganic_carbon_waste(i, j, k, grid,bgc.detritus, bgc, fields, auxiliary_fields)
    )
end

@inline function (bgc::OxygenNPD{<:Any, <:Nutrients{<:NitrateAmmonia}})(i, j, k, grid, ::Val{:O₂}, clock, fields, auxiliary_fields)
    rC = carbon_ratio(i, j, k, grid,bgc.plankton, bgc,fields)

    rP = bgc.oxygen.production_oxygen_carbon_ratio
    rN = bgc.oxygen.nitrification_oxygen_carbon_ratio

    # PP from NO₃ - nitrification etc.
    net_nitrate_production = bgc(i, j, k, grid, Val(:NO₃), clock, fields, auxiliary_fields)

    return (
        rP * primary_production(i, j, k, grid,bgc.plankton, bgc, fields, auxiliary_fields)
      - rN * rC * net_nitrate_production
      - rP * inorganic_carbon_waste(i, j, k, grid,bgc.plankton, bgc, fields, auxiliary_fields)
      - rP * inorganic_carbon_waste(i, j, k, grid,bgc.detritus, bgc, fields, auxiliary_fields)
    )
end

const OxygenCNDP{FT} = NutrientsPlanktonDetritus{FT, <:Nutrients{<:SingleTracerNutrient}, <:Any, <:CarbonNitrogenDissolvedParticulate, <:Any, <:Oxygen}

# we will have to think how to generalise this in the future for variable redfield planktons
@inline function (bgc::OxygenCNDP)(i, j, k, grid, ::Val{:O₂}, clock, fields, auxiliary_fields)
    rP = bgc.oxygen.production_oxygen_carbon_ratio
    rN = bgc.oxygen.nitrification_oxygen_carbon_ratio

    return (
        (rP + rN) * (primary_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
                     - inorganic_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields))
        - rP * inorganic_carbon_waste(i, j, k, grid,bgc.detritus, bgc, fields, auxiliary_fields)
        - rN * inorganic_nitrogen_waste(i, j, k, grid,bgc.detritus, bgc, fields, auxiliary_fields) * carbon_ratio(i, j, k, grid,bgc.plankton, bgc,fields) / nitrogen_ratio(i, j, k, grid,bgc.plankton, bgc,fields)
    )
end
end

module OxygenModels

export Oxygen

using ..NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus

using ..NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    NitrateAmmonia,
    nutrient_uptake,
    nitrification

using ..NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    primary_production,
    inorganic_carbon_waste

import Base: summary, show

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields

import ..NutrientsPlanktonDetritusModels: 
    carbon_ratio

"""
    Oxygen

`Oxygen` defines the evolution of oxygen (`O₂`) concentration for the
`LOBSTER` biogeochemical model. 

Oxygen is produced by photosynthesis in phytoplankton, and removed by nitrate 
production and oxidation of ammonia.

Oxygen concentration is only one way coupled with the rest of the biogeochemistry
and *does not* effect any other groups (e.g. low oxygen does *not* reduce
zooplankton growth). To capture this effect a different `plankton` could be defined.
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

# since we're assuming that nitrogen is nitrate (maybe erroniosuly), then it has to be denitrified before production/after remin
@inline function (bgc::OxygenNPD)(i, j, k, grid, ::Val{:O₂}, clock, fields, auxiliary_fields)
    rP = bgc.oxygen.production_oxygen_carbon_ratio
    rN = bgc.oxygen.nitrification_oxygen_carbon_ratio

    return (rP + rN) * (
        primary_production(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
      - inorganic_carbon_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
      - inorganic_carbon_waste(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
    )
end

@inline function (bgc::OxygenNPD{<:Any, <:Nutrients{<:NitrateAmmonia}})(i, j, k, grid, ::Val{:O₂}, clock, fields, auxiliary_fields)
    rC = carbon_ratio(bgc.plankton, bgc, i, j, k, fields)

    rP = bgc.oxygen.production_oxygen_carbon_ratio
    rN = bgc.oxygen.nitrification_oxygen_carbon_ratio

    # PP from NO₃ - nitrifcation etc.
    net_nitrate_production = bgc(i, j, k, grid, Val(:NO₃), clock, fields, auxiliary_fields)

    return (
        rP * primary_production(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
      - rN * rC * net_nitrate_production
      - rP * inorganic_carbon_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
      - rP * inorganic_carbon_waste(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
    )
end

end
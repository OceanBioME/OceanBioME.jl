"""
    Oxygen

`Oxygen` defines the evolution of oxygen (`O₂`) concentration for the
`LOBSTER` biogeochemical model. 

Oxygen is produced by photosynthesis in phytoplankton, and removed by nitrate 
production and oxidation of ammonia.

Oxygen concentration is only one way coupled with the rest of the biogeochemistry
and *does not* effect any other groups (e.g. low oxygen does *not* reduce
zooplankton growth). To capture this effect a different `biology` could be defined.
"""
@kwdef struct Oxygen{FT}
    respiration_oxygen_nitrogen_ratio :: FT = 10.75  # mol O/molN
  nitrification_oxygen_nitrogen_ratio :: FT = 2.0    # mol O/molN
end

required_biogeochemical_tracers(::Oxygen) = (:O₂, )

@inline function (lobster::LOBSTER{<:Any, <:Any, <:Any, <:Any, <:Oxygen})(i, j, k, grid, val_name::Val{:O₂}, clock, fields, auxiliary_fields)
    Rp = lobster.oxygen.respiration_oxygen_nitrogen_ratio
    Rn = lobster.oxygen.nitrification_oxygen_nitrogen_ratio

    μP = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    nitrate_production = lobster(i, j, k, grid, Val(:NH₄), clock, fields, auxiliary_fields)
    μNH₄ = nitrifcation(lobster.nutrients, i, j, k, fields) 

    return Rp * μP - (Rp - Rn) * nitrate_production - Rp * μNH₄
end

@inline function (lobster::LOBSTER{<:Any, <:Nutrient, <:Any, <:Any, <:Oxygen})(i, j, k, grid, val_name::Val{:O₂}, clock, fields, auxiliary_fields)
    Rp = lobster.oxygen.respiration_oxygen_nitrogen_ratio
    Rn = lobster.oxygen.nitrification_oxygen_nitrogen_ratio

    μP = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    respiration = biology_inorganic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)
                  + detritus_inorganic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)

    return Rp * μP - (Rp - Rn) * respiration
end

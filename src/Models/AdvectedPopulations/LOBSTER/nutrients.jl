"""
    NitrateAmmonia

`NitrateAmmonia` defines the default nutrient component for the `LOBSTER` 
biogeochemical model which includes nitrate (`NO₃`), and ammonia (`NH₄`)
in mmol N / m³.

Nitrate is only taken up by the biological component (by default only the
phytoplankton), and replenished by the nitrifcation of ammonia at the 
`nitrification_rate`.

Ammonia is also taken up by the biological component, and lost as it is 
turned into ammonia. It is replaced through waste from the biological system
(phytoplankton exudation, and zooplankton excretion), and from the degredation
of detritus.
"""
@kwdef struct NitrateAmmonia{FT}
    nitrification_rate::FT = 5.8e-7 # 1/s
end

required_biogeochemical_tracers(::NitrateAmmonia) = (:NO₃, :NH₄)

"""
    NitrateAmmoniaIron

`NitrateAmmoniaIron` defines a nutrient component for the `LOBSTER` 
biogeochemical model which includes nitrate (`NO₃`), and ammonia (`NH₄`)
in mmol N / m³ in the same way as the `NitrateAmmonia` component, but with
the addition of an iron component (`Fe`) in mmol Fe / m³.

Iron is only taken up by biological component and not replenished, as iron
gets remineralised a lot slower than the nitrogen and carbon in detritus so 
the iron is "never" returned to the biologically available pool. Instead it
is usually replenished by surface deposition of dust. For example, Signorini, 
et al. 2001 uses the same formulation.

Signorini, S. R., C. R. McClain, J. R. Christian, and C. S. Wong (2001), 
Seasonal and interannual variability of phytoplankton, nutrients, TCO2, pCO2, 
and O2 in the eastern subarctic Pacific (ocean weather station Papa), J. 
Geophys. Res., 106(C12), 31197–31215, doi:10.1029/2000JC000343.
"""
@kwdef struct NitrateAmmoniaIron{FT}
    nitrification_rate::FT = 5.8e-7 # 1/s
end

required_biogeochemical_tracers(::NitrateAmmoniaIron) = (:NO₃, :NH₄, :Fe)

@inline (lobster::LOBSTER{<:NitrateAmmoniaIron})(i, j, k, grid, val_name::Val{:Fe}, clock, fields, auxiliary_fields) =
    - nutrient_uptake(lobster, i, j, k, val_name, fields, auxiliary_fields)

##### common
const INCLUDES_NITRATE_AMMONIA = Union{NitrateAmmonia, NitrateAmmoniaIron}
const LOBSTER_WITH_NITRATE_AMMONIA = Union{LOBSTER{<:NitrateAmmonia}, LOBSTER{<:NitrateAmmoniaIron}}

@inline nitrifcation(nutrients::INCLUDES_NITRATE_AMMONIA, i, j, k, fields) = @inbounds nutrients.nitrification_rate * fields.NH₄[i, j, k]

@inline (lobster::LOBSTER_WITH_NITRATE_AMMONIA)(i, j, k, grid, val_name::Val{:NO₃}, clock, fields, auxiliary_fields) = (
    nitrifcation(lobster.nutrients, i, j, k, fields) 
  - nutrient_uptake(lobster, i, j, k, val_name, fields, auxiliary_fields)
) # done!

@inline (lobster::LOBSTER_WITH_NITRATE_AMMONIA)(i, j, k, grid, val_name::Val{:NH₄}, clock, fields, auxiliary_fields) = (
    biology_inorganic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)
  + detritus_inorganic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)
  - nitrifcation(lobster.nutrients, i, j, k, fields) 
  - nutrient_uptake(lobster, i, j, k, val_name, fields, auxiliary_fields)
) # done!
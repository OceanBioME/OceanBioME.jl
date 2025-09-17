"""
    CarbonateSystem

`CarbonateSystem` defines the carbonate system for the `LOBSTER` biogeochemical
model and evolves dissolved inorganic carbon (`DIC`) and alkalinity (`Alk`).
`DIC` is the total inorganic carbon, typically considered to be dissolved 
carbon dioxide, carbonic acid, bicarbonate, and carbonate. `Alk`alinity is 
the buffering capacity of the water which is approximatly the amount of bases
in solution that can be neutralised by the addition of acid so is approximatly
the sum of non-conserved ions such as carbonate, bicarbonate, hydroxide ions,
etc.

With this module activated the carbon budget can be closed as the flux of carbon
into the biology and detritus comes from the DIC. Additionally, the flux in/out 
of the water from the air can be computed via the `CarbonChemistry` module.

`DIC` is produced by remineralisation of organic matter, for example the breakdown
of detritus or excretion from zooplankton waste. It is also produced by the 
dissolution of calcite in the zooplankton gut. It is removed by photosynthesis
in phytoplankton.

`Alk`alinity is produced (by the removal of nitrose acid (?)), and removed by
the uptake of calcite into phytoplankton.

`DIC` and `Alk` concentration is only one way coupled with the rest of the 
biogeochemistry and *does not* effect any other groups (e.g. acidifcation does
not effect phytoplankton growth). To capture this effect a different `biology` 
could be defined.
"""
struct CarbonateSystem end

required_biogeochemical_tracers(::CarbonateSystem) = (:DIC, :Alk)

@inline (lobster::LOBSTER{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:DIC}, clock, fields, auxiliary_fields) = (
  - phytoplankton_primary_production(lobster, i, j, k, fields, auxiliary_fields)
  + biology_inorganic_carbon_waste(lobster, i, j, k, fields, auxiliary_fields)
  + detritus_inorganic_carbon_waste(lobster, i, j, k, fields, auxiliary_fields)
  + calcite_dissolution(lobster, i, j, k, fields, auxiliary_fields)
)

@inline (lobster::LOBSTER{<:Any, <:Any, <:Any, <:CarbonateSystem})(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields) = (
    nutrient_uptake(lobster, i, j, k, Val(:NO₃), fields, auxiliary_fields)
  - calcite_uptake(lobster, i, j, k, fields, auxiliary_fields)
)

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

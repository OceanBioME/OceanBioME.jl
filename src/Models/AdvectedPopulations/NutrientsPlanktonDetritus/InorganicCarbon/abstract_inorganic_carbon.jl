abstract type AbstractInorganicCarbon end

const NPD_AIC{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:Any, <:AbstractInorganicCarbon}

@inline (bgc::NPD_AIC)(i, j, k, grid, ::Val{:DIC}, clock, fields, auxiliary_fields) = (
  - primary_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_carbon_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
  - net_calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields)
)

@inline function (bgc::NPD_AIC)(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields)
    PN = phosphate_ratio(i, j, k, grid, bgc.plankton, bgc, fields) / nitrogen_ratio(i, j, k, grid, bgc.plankton, bgc, fields)

    return (
      - nitrate_primary_production(i, j, k, grid, bgc.plankton, bgc, clock, fields, auxiliary_fields) * (1 + PN)
      + ammonia_primary_production(i, j, k, grid, bgc.plankton, bgc, clock, fields, auxiliary_fields) * (1 - PN)
      - 2 * net_calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields)
    )
end

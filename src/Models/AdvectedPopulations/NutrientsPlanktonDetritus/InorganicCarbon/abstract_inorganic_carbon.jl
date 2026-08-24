abstract type AbstractInorganicCarbon end

const NPD_AIC{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:Any, <:AbstractInorganicCarbon}

@inline net_biological_dic_uptake(i, j, k, grid, bgc, fields, auxiliary_fields) = (
  - primary_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_carbon_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
)

@inline function net_biological_alkalinity_uptake(i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    PN = phosphate_ratio(i, j, k, grid, bgc.plankton, bgc, fields) / nitrogen_ratio(i, j, k, grid, bgc.plankton, bgc, fields)

    return (
      - nitrate_primary_production(i, j, k, grid, bgc.plankton, bgc, clock, fields, auxiliary_fields) * (1 + PN)
      + ammonia_primary_production(i, j, k, grid, bgc.plankton, bgc, clock, fields, auxiliary_fields) * (1 - PN)
    )
end

@inline (bgc::NPD_AIC)(i, j, k, grid, ::Val{:DIC}, clock, fields, auxiliary_fields) =
    net_biological_dic_uptake(i, j, k, grid, bgc, fields, auxiliary_fields) -
    net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields)

@inline (bgc::NPD_AIC)(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields) =
    net_biological_alkalinity_uptake(i, j, k, grid, bgc, clock, fields, auxiliary_fields) -
    2 * net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields)

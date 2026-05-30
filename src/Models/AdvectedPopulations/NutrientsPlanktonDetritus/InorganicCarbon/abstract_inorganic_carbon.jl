abstract type AbstractInorganicCarbon end

const NPD_AIC{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:Any, <:AbstractInorganicCarbon}

@inline (bgc::NPD_AIC)(i, j, k, grid, ::Val{:DIC}, clock, fields, auxiliary_fields) = (
  - primary_production(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
  + inorganic_carbon_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
  + inorganic_carbon_waste(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
  - net_calcite_production(bgc, i, j, k, fields, auxiliary_fields)
)

@inline function (bgc::NPD_AIC)(i, j, k, grid, ::Val{:Alk}, clock, fields, auxiliary_fields)
    PN = phosphate_ratio(bgc.plankton, bgc, i, j, k, fields) / nitrogen_ratio(bgc.plankton, bgc, i, j, k, fields)

    return (nitrate_primary_production(bgc.plankton, bgc, i, j, k, grid, clock, fields, auxiliary_fields) * (1 - PN)
            - ammonia_primary_production(bgc.plankton, bgc, i, j, k, grid, clock, fields, auxiliary_fields) * (1 + PN) 
            - 2 * net_calcite_production(bgc, i, j, k, fields, auxiliary_fields))
end

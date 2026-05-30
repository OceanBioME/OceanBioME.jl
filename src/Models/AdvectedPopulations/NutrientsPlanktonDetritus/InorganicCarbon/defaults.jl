@inline inorganic_carbon_waste(plankton_or_detritus, bgc, i, j, k,
                               fields, auxiliary_fields) =
    carbon_ratio(bgc.plankton, bgc, i, j, k, fields) * 
    inorganic_waste(plankton_or_detritus, bgc, i, j, k, fields, auxiliary_fields)

@inline dissolved_carbon_waste(plankton_or_detritus, bgc, i, j, k,
                               fields, auxiliary_fields) =
    carbon_ratio(bgc.plankton, bgc, i, j, k, fields) * 
    dissolved_waste(plankton_or_detritus, bgc, i, j, k, fields, auxiliary_fields)

@inline primary_production(plankton, bgc, i, j, k, fields, auxiliary_fields) =
    carbon_ratio(bgc.plankton, bgc, i, j, k, fields) * 
    nutrient_uptake(plankton, bgc, i, j, k, fields, auxiliary_fields)

# assume PP is nitrate based if nitrate and ammonia not available
@inline nitrate_primary_production(plankton,
                                   bgc::NPD{FT, <:Nutrients{SingleTracerNutrient}}, 
                                   i, j, k, grid, clock, fields, auxiliary_fields) where FT = 
    bgc(i, j, k, grid, Val(:N), clock, fields, auxiliary_fields)

@inline nitrate_primary_production(plankton,
                                   bgc::NPD{FT, <:Nutrients{<:NitrateAmmonia}}, 
                                   i, j, k, grid, clock, fields, auxiliary_fields) where FT = 
    bgc(i, j, k, grid, Val(:NO₃), clock, fields, auxiliary_fields)

@inline ammonia_primary_production(plankton, ::NPD{FT}, args...) where FT = zero(FT)

@inline ammonia_primary_production(plankton,
                                   bgc::NPD{FT, <:Nutrients{<:NitrateAmmonia}}, 
                                   i, j, k, grid, clock, fields, auxiliary_fields) where FT = 
    bgc(i, j, k, grid, Val(:NH₄), clock, fields, auxiliary_fields)

# default assumption that the biological and detrital pool stores calcite at a fixed ratio to carbon
# but the dissolved elements might not so `calcite_dissolution` might not be = R(1+ρ)*inorganic_waste
@inline net_calcite_production(bgc, i, j, k, fields, auxiliary_fields) = (
    calcite_rain_ratio(bgc.plankton, bgc, i, j, k, fields)  * (
        primary_production(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
      - inorganic_carbon_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
    ) 
  - calcite_dissolution(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
)

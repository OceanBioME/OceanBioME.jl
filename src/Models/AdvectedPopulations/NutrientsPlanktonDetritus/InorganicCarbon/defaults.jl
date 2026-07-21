@inline primary_production(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
    carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields) * 
    nutrient_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

# assume PP is nitrate based if nitrate and ammonia not available
@inline nitrate_primary_production(i, j, k, grid, 
                                   plankton, bgc::NPD{FT, <:Nutrients{SingleTracerNutrient}}, 
                                   clock, fields, auxiliary_fields) where FT = 
    bgc(i, j, k, grid, Val(:N), clock, fields, auxiliary_fields)

@inline nitrate_primary_production(i, j, k, grid, 
                                   plankton, bgc::NPD{FT, <:Nutrients{Nothing}}, 
                                   clock, fields, auxiliary_fields) where FT = 
    nitrogen_ratio(i, j, k, grid, plankton, bgc, fields) * (
            inorganic_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
          + inorganic_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
          - nutrient_uptake(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
    )

@inline nitrate_primary_production(i, j, k, grid, 
                                   plankton, bgc::NPD{FT, <:Nutrients{<:NitrateAmmonia}}, 
                                   clock, fields, auxiliary_fields) where FT = 
    bgc(i, j, k, grid, Val(:NO₃), clock, fields, auxiliary_fields)

@inline ammonia_primary_production(i, j, k, grid, plankton, ::NPD{FT}, args...) where FT = zero(FT)

@inline ammonia_primary_production(i, j, k, grid, 
                                   plankton, bgc::NPD{FT, <:Nutrients{<:NitrateAmmonia}}, 
                                   clock, fields, auxiliary_fields) where FT = 
    bgc(i, j, k, grid, Val(:NH₄), clock, fields, auxiliary_fields)

# default assumption that the biological and detrital pool stores calcite at a fixed ratio to carbon
# but the dissolved elements might not so `calcite_dissolution` might not be = R(1+ρ)*inorganic_waste
@inline net_calcite_production(i, j, k, grid, bgc, fields, auxiliary_fields) = (
    calcite_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)  * (
        primary_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
      - inorganic_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
    ) 
  - calcite_dissolution(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
)

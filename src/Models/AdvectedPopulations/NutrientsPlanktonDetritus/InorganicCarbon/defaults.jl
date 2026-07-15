@inline primary_production(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
    carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields) *
    nutrient_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

# Default methods for the three per-plankton calcium carbonate hooks owned by this module (used by
# ExplicitCalciumCarbonate; PhytoZoo overrides all three in Plankton/phyto_zoo.jl). They hold for any plankton
# with no standing biomass tracer and together conserve carbon exactly there, because primary
# production then equals the sum of the dissolved + inorganic + solid carbon waste:
#   - biological_calcium_carbonate_precipitation — DIC-side formation: ρ × primary production
#   - particulate_calcium_carbonate_production   — into the solid CaCO₃ pool: ρ × solid carbon waste
#   - biological_calcium_carbonate_dissolution  — returned straight to DIC/Alk: ρ × (inorganic + dissolved) waste
@inline biological_calcium_carbonate_precipitation(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
    calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields) *
    primary_production(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

@inline particulate_calcium_carbonate_production(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
    calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields) *
    solid_carbon_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

@inline biological_calcium_carbonate_dissolution(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
    calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields) * (
        inorganic_carbon_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
      + dissolved_carbon_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
    )

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

@inline net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields) = (
    calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)  * (
        primary_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
      - inorganic_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
    ) 
  - calcium_carbonate_dissolution(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
)

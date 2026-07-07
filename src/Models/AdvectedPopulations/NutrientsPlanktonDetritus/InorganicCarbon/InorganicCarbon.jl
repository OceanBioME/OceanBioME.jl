module InorganicCarbonModels

export CarbonateSystem, ExplicitCalcite

using Oceananigans.Units

using ..NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus,
    NPD,
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio,
    calcite_rain_ratio

using ..NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    SingleTracerNutrient,
    NitrateAmmonia

using ..NutrientsPlanktonDetritusModels:
    dissolved_waste,
    calcite_dissolution,
    nutrient_uptake,
    inorganic_waste,
    inorganic_carbon_waste,
    solid_carbon_waste,
    dissolved_carbon_waste

import Base: summary, show

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields

import ..NutrientsPlanktonDetritusModels:
    carbon_ratio

# the per-plankton calcite hooks (biological_calcite_precipitation, particulate_calcite_production,
# biological_calcite_dissolution) are owned by this module — their default methods in defaults.jl
# create the generic functions; PhytoZoo overrides them in Plankton/phyto_zoo.jl.

include("abstract_inorganic_carbon.jl")
include("defaults.jl")
include("implicit_calcite.jl")
include("explicit_calcite.jl")

end # module
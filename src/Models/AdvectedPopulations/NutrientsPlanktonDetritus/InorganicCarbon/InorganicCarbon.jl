module InorganicCarbonModels

export CarbonateSystem, ExplicitCalciumCarbonate

using Oceananigans.Units

using ..NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus,
    NPD,
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio,
    calcium_carbonate_rain_ratio

using ..NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    SingleTracerNutrient,
    NitrateAmmonia

using ..NutrientsPlanktonDetritusModels:
    dissolved_waste,
    calcium_carbonate_dissolution,
    nutrient_uptake,
    nitrogen_fixation,
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

include("abstract_inorganic_carbon.jl")
include("defaults.jl")
include("implicit_calcium_carbonate.jl")
include("explicit_calcium_carbonate.jl")

end # module
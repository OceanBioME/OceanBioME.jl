module InorganicCarbonModels

export CarbonateSystem, ExplicitCalciumCarbonate, ImplicitExplicitCalcite

using Adapt
using Oceananigans.Units
using Oceananigans.Grids: AbstractGrid, Center
using Oceananigans.Fields: CenterField, Field

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
    inorganic_waste,
    inorganic_carbon_waste,
    solid_carbon_waste,
    dissolved_carbon_waste

import Adapt: adapt_structure, adapt
import Base: summary, show

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_auxiliary_fields,
    update_biogeochemical_state!

import ..NutrientsPlanktonDetritusModels:
    carbon_ratio,
    group_element_tracers

include("abstract_inorganic_carbon.jl")
include("defaults.jl")
include("implicit_calcium_carbonate.jl")
include("explicit_calcium_carbonate.jl")
include("implicit_explicit_calcite.jl")

end # module
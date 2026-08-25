module PlanktonModels

export Abiotic, ImplicitProductivity, PhytoZoo

# plankton models *must* define `inorganic_waste`, `nutrient_uptake`, `dissolved_waste`, and `solid_waste`
# and may define `carbon_ratio`, `nitrogen_ratio`, `phosphate_ratio`, `iron_ratio`, `silicon_ratio`, and `calcite_rain_ratio`
# if they are not defined they use the default elemental ratios, with nitrogen as the base unit
# `chlorophyll_ratio` may be defined once for the component and overridden for individual tracers
# you may also define `X_Y_waste` where X is `inorganic`, `dissolved` and `solid`, and `Y` are the elements

# plankton may also define `grazing` for detritus tracers they consume; it otherwise defaults to zero
# and can represent zooplankton grazing, dissolved uptake like MARBL, or mixotrophy

# `nutrient_uptake` is called with `i, j, k, grid, plankton, bgc, fields, auxiliary_fields`,
# and may also be defined with `Val(name)` between `grid` and `plankton` for tracer-specific uptake

using Adapt
using Oceananigans.Units
using Oceananigans.Grids: AbstractGrid
using OceanBioME: setup_velocity_fields

using ..NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus,
    NPD

using ..NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    SingleTracerNutrient,
    NitrateAmmonia

using ..NutrientsPlanktonDetritusModels.DetritusModels:
    Detritus,
    DissolvedParticulate,
    InstantRemineralisationDetritus,
    CarbonNitrogenDissolvedParticulate

import Adapt: adapt_structure
import Base: show, summary

import Oceananigans.Biogeochemistry: 
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields

import OceanBioME: chlorophyll

import ..NutrientsPlanktonDetritusModels:
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio,
    calcium_carbonate_rain_ratio

import ..NutrientsPlanktonDetritusModels:
    inorganic_waste,
    nutrient_uptake,
    solid_waste,
    dissolved_waste,
    inorganic_nitrogen_waste,
    inorganic_phosphate_waste,
    inorganic_iron_waste,
    inorganic_silicon_waste,
    solid_nitrogen_waste,
    solid_phosphate_waste,
    solid_iron_waste,
    solid_silicon_waste,
    dissolved_nitrogen_waste,
    dissolved_phosphate_waste,
    dissolved_iron_waste,
    dissolved_silicon_waste

import ..NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    inorganic_carbon_waste,
    primary_production,
    net_calcium_carbonate_production,
    calcium_carbonate_rain_ratio,
    biological_calcium_carbonate_precipitation,
    particulate_calcium_carbonate_production,
    biological_calcium_carbonate_dissolution

import ..NutrientsPlanktonDetritusModels.DetritusModels:
    grazing,
    calcium_carbonate_precipitation

chlorophyll(bgc::NutrientsPlanktonDetritus, model) =
    chlorophyll(bgc.plankton, model)

include("primitives.jl")
include("abiotic.jl")
include("implicit.jl")
include("phyto_zoo.jl")

end # module
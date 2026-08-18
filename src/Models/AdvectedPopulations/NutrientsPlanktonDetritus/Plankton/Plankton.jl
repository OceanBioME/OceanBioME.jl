module PlanktonModels

export Abiotic, ImplicitProductivity, PhytoZoo

# Plankton components define `inorganic_waste`, `nutrient_uptake`, `dissolved_waste`, and `solid_waste`.
# Elemental composition is provided by `carbon_ratio`, `nitrogen_ratio`, `phosphate_ratio`, `iron_ratio`,
# `silicon_ratio`, and `calcite_rain_ratio`; the NPD defaults use nitrogen as the base unit.
# Components may specialise `plankton_element_tracers` when their owned tracers require conservation
# coefficients that differ from the community-wide elemental ratios.
#
# Detritus consumption uses the generic
# `grazing(i, j, k, grid, Val(detritus_name), plankton, bgc, fields, auxiliary_fields)` dispatch.
# The default is zero, while plankton components may specialise the method for any realised detritus name.
#
# `nutrient_uptake` may be defined as an aggregate uptake or specialised for `Val(nutrient_name)`.

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
    calcite_rain_ratio

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
    net_calcite_production,
    calcite_rain_ratio

import ..NutrientsPlanktonDetritusModels.DetritusModels:
    grazing,
    calcite_precipitation

chlorophyll(bgc::NutrientsPlanktonDetritus, model) =
    chlorophyll(bgc.plankton, model)

include("primitives.jl")
include("abiotic.jl")
include("implicit.jl")
include("phyto_zoo.jl")

end # module
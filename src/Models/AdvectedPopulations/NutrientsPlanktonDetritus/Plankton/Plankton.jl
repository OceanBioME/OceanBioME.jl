module PlanktonModels

export Abiotic, ImplicitProductivity, PhytoZoo

# Plankton components define `inorganic_waste`, `nutrient_uptake`, `dissolved_waste`, and `solid_waste`.
# Elemental composition is provided by `carbon_ratio`, `nitrogen_ratio`, `phosphate_ratio`, `iron_ratio`,
# `silicon_ratio`, and `calcite_rain_ratio`; the defaults use nitrogen as the base unit.
# `plankton_element_tracers` can provide tracer-specific composition for conservation calculations.
#
# Plankton components can consume detritus by defining `grazing` for the detritus tracers they consume.
# `nutrient_uptake` may be defined for the whole plankton component or for an individual nutrient tracer.

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
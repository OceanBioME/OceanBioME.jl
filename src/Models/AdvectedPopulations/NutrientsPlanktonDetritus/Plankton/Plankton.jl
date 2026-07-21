module PlanktonModels

export Abiotic, ImplicitProductivity, PhytoZoo

# plankton models *must* define `inorganic_waste`, `nutrient_uptake`, `dissolved_waste`, and `solid_waste`
# and may define `carbon_ratio`, `nitrogen_ratio`, `phosphate_ratio`, `iron_ratio`, and `silicon_ratio`
# if they are not defined they default to 106/16:1:0.0062/16:0 (i.e. the default unit is nitrogen)
# you may also define `X_Y_waste` where X is `inorganic`, `dissolved` and `solid`, and `Y` are the elements

# plankton may also define `detritus_grazing(plankton, bgc, i, j, k, val_name, fields, auxiliary_fields)` which
# acts on detritus pools and otherwise defaults to zero (and is representative of both zooplankton grazing,
# dissolved uptake like MARBL, or mixotrophy)

# they must define `nutrient_uptake` with arguments `plankton, bgc, i, j, k, fields, auxiliary_fields`,
# but may specialise to arguments `plankton, bgc, i, j, k, val_name, fields, auxiliary_fields`

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
module DetritusModels

export Detritus, DissolvedParticulate, InstantRemineralisationDetritus, CarbonNitrogenDissolvedParticulate

using Adapt
using Oceananigans.Grids: AbstractGrid
using OceanBioME: setup_velocity_fields

using ..NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus,
    NPD,
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio,
    calcite_rain_ratio

import ..NutrientsPlanktonDetritusModels: dissolved_waste, solid_waste, calcite_dissolution, inorganic_waste, nutrient_uptake

using ..NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    SingleTracerNutrient,
    NitrateAmmonia

import Adapt: adapt_structure
import Base: summary, show

import Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity

import ..NutrientsPlanktonDetritusModels:
    inorganic_nitrogen_waste,
    inorganic_phosphate_waste,
    inorganic_iron_waste,
    inorganic_silicon_waste

include("defaults.jl")
include("instant_remineralisation.jl")
include("single_detritus.jl")
include("single_element.jl")
include("carbon_nitrogen.jl")

end # module
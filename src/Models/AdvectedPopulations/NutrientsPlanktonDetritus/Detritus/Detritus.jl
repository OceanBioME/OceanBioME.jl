module DetritusModels

export Detritus, DissolvedParticulate, InstantRemineralisation

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

import ..NutrientsPlanktonDetritusModels.NutrientsModels: 
    inorganic_waste,
    inorganic_nitrogen_waste,
    inorganic_phosphate_waste,
    inorganic_iron_waste,
    inorganic_silicon_waste,
    nutrient_uptake

include("defaults.jl")
include("instant_remineralisation.jl")
include("single_detritus.jl")
include("single_element.jl")

end # module
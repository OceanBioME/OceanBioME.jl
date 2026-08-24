module MultiElement

export MultiElementRefractoryDissolvedParticulate, 
       MultiElementRefractoryDissolved

using Adapt
using Oceananigans.Units: day
using Oceananigans.Grids: AbstractGrid, znode, Center, Face
using Oceananigans.Operators: Δzᶜᶜᶜ
using Oceananigans.Fields: CenterField, ZFaceField, Field
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell
using Oceananigans.Architectures: architecture
using Oceananigans.Models: fields
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index
using OceanBioME: setup_velocity_fields

using ...NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus,
    NPD,
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio,
    calcium_carbonate_rain_ratio

import ...NutrientsPlanktonDetritusModels: 
    dissolved_waste, solid_waste, 
    calcium_carbonate_dissolution, 
    inorganic_waste, nutrient_uptake

using ...NutrientsPlanktonDetritusModels.NutrientsModels:
    Nutrients,
    SingleTracerNutrient,
    NitrateAmmonia

import Adapt: adapt_structure
import Base: summary, show

import Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity,
    update_biogeochemical_state!

import ...NutrientsPlanktonDetritusModels:
    inorganic_nitrogen_waste,
    inorganic_phosphate_waste,
    inorganic_iron_waste,
    inorganic_silicon_waste

using ...NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    particulate_calcium_carbonate_production

using OceanBioME.Light: @subcolumn_average, @preserve_subcolumns, interface_par

using ...NutrientsPlanktonDetritusModels.NutrientsModels:
    ComplexedIron,
    ballast_sinking_mass

import ...NutrientsPlanktonDetritusModels.NutrientsModels:
    poc_remineralisation,
    doc_production,
    scavenging_sinking_mass

function dissolved_phosphorus_production end
function particulate_phosphorus_production end
function dissolved_phosphorus_balance end
function dissolved_phosphorus_uptake end
function particulate_iron_production end
function silicon_dissolution end
function biogenic_silica_production end

function ballast_oxygen_scale end

const yr = 365day

include("ballast.jl")
include("explicit_particles.jl")
include("implicit_particles.jl")
include("adapt_show_methods.jl")

end # module
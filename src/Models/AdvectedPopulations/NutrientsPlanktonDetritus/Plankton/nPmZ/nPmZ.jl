# MARBLS nPmZ plankton model where phytoplankton pools have carbon,
# chlorophyll, phosphorus, iron, and optionally silicon and calcite
# and zooplankton pools have carbon only
module MARBLPlankton # maybe rename

export ManyPhytoZoo, Autotroph, Heterotroph, Grazing
export MARBL_small_phyto, MARBL_diatoms, MARBL_diazotrophs, MARBL_autotrophs, MARBL_zooplankton
export MARBL_cocco_small_phyto, MARBL_cocco_diatoms, MARBL_cocco_diazotrophs, MARBL_coccolithophores
export MARBL_cocco_autotrophs, MARBL_cocco_zooplankton, MARBL_cocco_plankton

# owned here and used by the components that couple to the plankton (see the second include pass in
# NutrientsPlanktonDetritus.jl)
export oxygen_production_total

using Adapt
using Oceananigans: defaults
using Oceananigans.Units: day
using Oceananigans.Fields: ZeroField, CenterField, Field, AbstractField
using Oceananigans.Grids: znode, Center, Face, AbstractGrid
using Oceananigans.Operators: ℑzᵃᵃᶠ, Δzᶜᶜᶜ
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Models: fields
using Oceananigans.Utils: launch!
using KernelAbstractions: @kernel, @index

# the enclosing modules are still being defined at include time, so every reference
# below is relative (`..` = PlanktonModels, `...` = NutrientsPlanktonDetritusModels)
using ..PlanktonModels: AbstractPlankton

using ...NutrientsPlanktonDetritusModels: NutrientsPlanktonDetritus, NPD
using ...NutrientsPlanktonDetritusModels.NutrientsModels: Nutrients, NitrateAmmonia
using ...NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    AbstractInorganicCarbon,
    ExplicitCalciumCarbonate,
    carbonate_field_names,
    hydrostatic_pressure

using OceanBioME.Light: @subcolumn_average, @preserve_subcolumns

import ...NutrientsPlanktonDetritusModels.NutrientsModels:
    iron_scavenging,
    scavenging_sinking_mass

import ...NutrientsPlanktonDetritusModels.DetritusModels.MultiElement:
    dissolved_phosphorus_production,
    particulate_phosphorus_production,
    dissolved_phosphorus_balance,
    dissolved_phosphorus_uptake,
    particulate_iron_production,
    silicon_dissolution,
    biogenic_silica_production

import Adapt: adapt_structure
import Base: summary, show

import Oceananigans.Biogeochemistry:
    required_biogeochemical_tracers,
    required_biogeochemical_auxiliary_fields,
    biogeochemical_drift_velocity,
    update_biogeochemical_state!

import OceanBioME: chlorophyll, conserved_tracers

import ...NutrientsPlanktonDetritusModels:
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio,
    nutrient_uptake,
    inorganic_waste,
    solid_waste,
    dissolved_waste,
    inorganic_carbon_waste,
    inorganic_nitrogen_waste,
    inorganic_phosphate_waste,
    inorganic_iron_waste,
    group_element_tracers

# detrital calcite is rehomed onto ExplicitCalciumCarbonate: we override its three per-plankton
# calcite hooks (summed over the calcifying PFTs) rather than net_calcium_carbonate_production
import ...NutrientsPlanktonDetritusModels.InorganicCarbonModels:
    primary_production,
    biological_calcium_carbonate_precipitation,
    particulate_calcium_carbonate_production,
    biological_calcium_carbonate_dissolution,
    calcium_carbonate_dissolution_flux

include("autotroph.jl")
include("heterotroph.jl")
include("primatives.jl")
include("many_phyto_zoo.jl")
include("manifest_many_phyto_zoo.jl")
include("infrastructure.jl")
include("grazing.jl")
include("marbl_pfts.jl")
include("defaults.jl")
include("adapt_show_methods.jl")

end # module

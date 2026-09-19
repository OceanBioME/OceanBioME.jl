module DetritusModels

export Detritus, DissolvedParticulate, InstantRemineralisationDetritus, CarbonNitrogenDissolvedParticulate

using Adapt
using Oceananigans.Grids: AbstractGrid
using OceanBioME: setup_velocity_fields, ExplicitSinking, ImplicitSinking, floor_index_field,
                  implicit_sinking_column!, dissolution_length

import OceanBioME: implicit_sinking_production

using ..NutrientsPlanktonDetritusModels:
    NutrientsPlanktonDetritus,
    NPD,
    carbon_ratio,
    nitrogen_ratio,
    phosphate_ratio,
    iron_ratio,
    silicon_ratio,
    calcium_carbonate_rain_ratio

import ..NutrientsPlanktonDetritusModels: dissolved_waste, solid_waste, calcium_carbonate_dissolution, inorganic_waste, nutrient_uptake

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
    biogeochemical_drift_velocity,
    update_biogeochemical_state!

using Oceananigans.Architectures: architecture
using Oceananigans.Utils: launch!
using Oceananigans: fields

import ..NutrientsPlanktonDetritusModels:
    inorganic_nitrogen_waste,
    inorganic_phosphate_waste,
    inorganic_iron_waste,
    inorganic_silicon_waste

"""
    AbstractSinkingDetritus{SK}

Abstract supertype for detritus models that carry a `sinking :: SK` field, where `SK` is either
[`ExplicitSinking`](@ref) or [`ImplicitSinking`](@ref). Provides a shared
`update_biogeochemical_state!` implementation for implicit sinking.
"""
abstract type AbstractSinkingDetritus{SK} end

include("defaults.jl")
include("instant_remineralisation.jl")
include("single_detritus.jl")
include("single_element.jl")
include("carbon_nitrogen.jl")

# --- shared implicit sinking update for all AbstractSinkingDetritus ---

function update_biogeochemical_state!(model, detritus::AbstractSinkingDetritus{<:ImplicitSinking}, npd::NutrientsPlanktonDetritus)
    sinking = detritus.sinking
    grid = model.grid
    Nz = size(grid, 3)
    FT = eltype(grid)

    for name in keys(sinking.remineralisation)
        ℓ = dissolution_length(sinking, name)
        launch!(architecture(grid), grid, :xy, implicit_sinking_column!,
                grid, detritus, npd, fields(model), biogeochemical_auxiliary_fields(model.biogeochemistry),
                sinking.remineralisation[name], sinking.floor_flux[name],
                sinking.floor_indices, convert(FT, ℓ),
                sinking.open_bottom, Nz, Val(name))
    end

    return nothing
end

end # module

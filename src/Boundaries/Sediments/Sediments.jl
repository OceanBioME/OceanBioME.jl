module Sediments

export SimpleMultiG, InstantRemineralisation

using KernelAbstractions

using OceanBioME: Biogeochemistry, BoxModelGrid

using Oceananigans
using Oceananigans.Architectures: device, architecture, arch_array
using Oceananigans.Utils: launch!
using Oceananigans.Advection: advective_tracer_flux_z
using Oceananigans.Units: day
using Oceananigans.Fields: ConstantField
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Grids: zspacing
using Oceananigans.Operators: volume
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell


import Adapt: adapt_structure, adapt
import Oceananigans.Biogeochemistry: update_tendencies!

abstract type AbstractSediment end
abstract type FlatSediment <: AbstractSediment end

@inline bottom_index(i, j, sediment) = @inbounds bottom_index_array(sediment)[i, j]

sediment_fields(::AbstractSediment) = ()

function update_tendencies!(bgc, sediment::FlatSediment, model)
    arch = model.grid.architecture

    for (i, tracer) in enumerate(sediment_tracers(sediment))
        launch!(arch, model.grid, :xy, store_flat_tendencies!, sediment.tendencies.G⁻[i], sediment.tendencies.Gⁿ[i])
    end

    launch!(arch, model.grid, :xy,
            _calculate_tendencies!,
            sediment, bgc.underlying_biogeochemistry, model.grid, sinking_advection(bgc.underlying_biogeochemistry, model.advection), model.tracers, model.timestepper.Gⁿ, sediment.tendencies.Gⁿ)
            
    return nothing
end


@kernel function store_flat_tendencies!(G⁻, G⁰)
    i, j = @index(Global, NTuple)
    @inbounds G⁻[i, j, 1] = G⁰[i, j, 1]
end

@inline nitrogen_flux() = 0
@inline carbon_flux() = 0
@inline remineralisation_receiver() = nothing
@inline sinking_tracers() = nothing

@inline sinking_advection(bgc, advection) = advection
@inline sinking_advection(bgc, advection::NamedTuple) = advection[sinking_tracers(bgc)]

@inline nitrogen_flux(i, j, k, grid, advection, bgc::Biogeochemistry, tracers) = nitrogen_flux(i, j, k, grid, advection, bgc.underlying_biogeochemistry, tracers)
@inline carbon_flux(i, j, k, grid, advection, bgc::Biogeochemistry, tracers) = carbon_flux(i, j, k, grid, advection, bgc.underlying_biogeochemistry, tracers)
@inline remineralisation_receiver(bgc::Biogeochemistry) = remineralisation_receiver(bgc.underlying_biogeochemistry)

@inline advection_scheme(advection, val_tracer) = advection
@inline advection_scheme(advection::NamedTuple, val_tracer::Val{T}) where T = advection[T]

@inline function sinking_flux(i, j, k, grid, advection, val_tracer::Val{T}, bgc, tracers) where T 
    return - advective_tracer_flux_z(i, j, k, grid, advection_scheme(advection, val_tracer), biogeochemical_drift_velocity(bgc, val_tracer).w, tracers[T]) /
      volume(i, j, k, grid, Center(), Center(), Center())
end

calculate_bottom_indices(::BoxModelGrid) = 1
calculate_bottom_indices(grid) = ones(Int, size(grid)[1:2]...)

@kernel function find_bottom_cell(grid, bottom_indices)
    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)

    k_bottom = 1

    while immersed_cell(i, j, k_bottom, grid.underlying_grid, grid.immersed_boundary) && k_bottom < Nz
        k_bottom += 1
    end

    @inbounds bottom_indices[i, j] = k_bottom
end

function calculate_bottom_indices(grid::ImmersedBoundaryGrid)
    arch = architecture(grid)
    indices = arch_array(arch, zeros(Int, size(grid)[1:2]...))

    launch!(grid.architecture, grid, :xy, find_bottom_cell, grid, indices)

    return indices
end

include("coupled_timesteppers.jl")
include("simple_multi_G.jl")
include("instant_remineralization.jl")

end # module

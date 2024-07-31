using Adapt

using Base: @propagate_inbounds
using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Grids: znode
using Oceananigans.Utils: launch!

import Adapt: adapt_structure
import Base: getindex
import Oceananigans.Biogeochemistry: update_biogeochemical_state!

# common functions for all mixed layer depth computation 
abstract type AbstractMixedLayerDepth end

@propagate_inbounds Base.getindex(f::AbstractMixedLayerDepth, inds...) = getindex(f.field, inds...)

function update_biogeochemical_state!(model, mld::AbstractMixedLayerDepth)
    grid = model.grid

    arch = architecture(grid)

    launch!(arch, grid, :xy, _compute_mixed_layer_depth!, mld, grid, model.tracers[required_tracers(mld)]...)

    fill_halo_regions!(mld.field)

    return nothing
end

@kwdef struct TemperatureChangeMixedLayerDepth{F, FT} <: AbstractMixedLayerDepth
                        field :: F
  temperature_change_criteria :: FT = 0.8 # Kara et. al. (2000)
end

TemperatureChangeMixedLayerDepth(; grid, temperature_change_criteria = 0.125) =
    TemperatureChangeMixedLayerDepth(Field{Center, Center, Nothing}(grid; indices = (:, :, 1:1)), 
                                    temperature_change_criteria)

Adapt.adapt_structure(to, mld::TemperatureChangeMixedLayerDepth) =
    TemperatureChangeMixedLayerDepth(adapt(to, mld.field),
                                     adapt(to, mld.temperature_change_criteria))

required_tracers(::TemperatureChangeMixedLayerDepth) = (:T, )

@kernel function _compute_mixed_layer_depth!(mld::TemperatureChangeMixedLayerDepth, grid, T)
    i, j = @index(Global, NTuple)

    k_surface = grid.Nz

    k = 1
    found = false
    for k′ in 1:k_surface # maybe this is unrollable 
        k = ifelse(T[i, j, k′] < T[i, j, k_surface] - mld.temperature_change_criteria, k′, k)
        found = ifelse(!found, T[i, j, k′] < T[i, j, k_surface] - mld.temperature_change_criteria, found)
    end

    @inbounds mld.field[i, j] = ifelse(found, znode(i, j, k, grid, Center(), Center(), Center()), -Inf)
end
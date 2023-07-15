module Sediments

export SimpleMultiG, InstantRemineralisation

using KernelAbstractions
using OceanBioME: ContinuousFormBiogeochemistry
using Oceananigans
using Oceananigans.Architectures: device
using Oceananigans.Utils: launch!
using Oceananigans.Advection: advective_tracer_flux_z
using Oceananigans.Units: day
using Oceananigans.Fields: CenterField, Face
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Grids: zspacing
using Oceananigans.Operators: volume
using Oceananigans.Fields: Center

import Adapt: adapt_structure, adapt

abstract type AbstractSediment end
abstract type FlatSediment <: AbstractSediment end

sediment_fields(::AbstractSediment) = ()

@inline update_sediment_tendencies!(bgc, sediment, model) = nothing

function update_sediment_tendencies!(bgc, sediment::FlatSediment, model)
    arch = model.grid.architecture

    for (i, tracer) in enumerate(sediment_tracers(sediment))
        launch!(arch, model.grid, :xy, store_flat_tendencies!, sediment.tendencies.G⁻[i], sediment.tendencies.Gⁿ[i])
    end

    launch!(arch, model.grid, :xy,
            _calculate_tendencies!,
            bgc.sediment_model, bgc, model.grid, model.advection, model.tracers, model.timestepper)
            
    return nothing
end

@kernel function store_flat_tendencies!(G⁻, G⁰)
    i, j = @index(Global, NTuple)
    @inbounds G⁻[i, j, 1] = G⁰[i, j, 1]
end


@inline nitrogen_flux(grid, adveciton, bgc, tracers, i, j) = 0
@inline carbon_flux(grid, adveciton, bgc, tracers, i, j) = 0
@inline remineralisation_receiver(bgc, tendencies) = nothing

@inline sinking_flux(i, j, grid, advection, val_tracer::Val{T}, bgc, tracers) where T = 
    - advective_tracer_flux_z(i, j, 1, grid, advection, biogeochemical_drift_velocity(bgc, val_tracer).w, tracers[T]) /
      volume(i, j, 1, grid, Center(), Center(), Center())

include("coupled_timesteppers.jl")
include("simple_multi_G.jl")
include("instant_remineralization.jl")

end # module

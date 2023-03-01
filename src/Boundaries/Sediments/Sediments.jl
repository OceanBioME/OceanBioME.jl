module Sediments

export SimpleMultiG

using KernelAbstractions
using OceanBioME: ContinuousFormBiogeochemistry
using Oceananigans.Architectures: device_event, device
using Oceananigans.Utils: launch!
using Oceananigans.Advection: div_Uc
using Oceananigans.Units: day
using Oceananigans.Fields: CenterField, Face
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity, biogeochemical_advection_scheme
using Oceananigans.Grids: znode

import Oceananigans.Biogeochemistry: update_tendencies!

abstract type AbstractSediment end
abstract type FlatSediment <: AbstractSediment end

sediment_fields(::AbstractSediment) = ()

function update_tendencies!(bgc::ContinuousFormBiogeochemistry{<:Any, <:FlatSediment}, model)
    arch = model.grid.architecture

    events = []

    sediment = bgc.sediment_model

    for (i, tracer) in enumerate(sediment_tracers(sediment))    
        field_event = launch!(arch, model.grid, :xy, store_flat_tendencies!, sediment.tendencies.Gⁿ[i], sediment.tendencies.G⁻[i], dependencies = device_event(arch))

        push!(events, field_event)
    end

    wait(device(model.architecture), MultiEvent(Tuple(events)))

    event = launch!(arch, model.grid, :xy,
                    calculate_tendencies!,
                    bgc.sediment_model, bgc, model,
                    dependencies = device_event(arch))

    wait(device(arch), event)
    return nothing
end

@kernel function store_flat_tendencies!(G⁻, G⁰)
    i, j = @index(Global, NTuple)
    @inbounds G⁻[i, j, 1] = G⁰[i, j, 1]
end

include("coupled_timesteppers.jl")
include("simple_multi_G.jl")

end # module
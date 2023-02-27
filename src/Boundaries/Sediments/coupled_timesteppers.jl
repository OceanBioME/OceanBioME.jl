using Oceananigans: NonhydrostaticModel, prognostic_fields
using OceanBioME: ContinuousFormBiogeochemistry
using OceanBioME.Boundaries.Sediments: AbstractSediment
using Oceananigans.TimeSteppers: ab2_step_field!
using Oceananigans.Utils: work_layout, launch!
using Oceananigans.Architectures: device_event
using Oceananigans.TurbulenceClosures: implicit_step!

import Oceananigans.TimeSteppers: ab2_step!

function ab2_step!(model::NonhydrostaticModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:ContinuousFormBiogeochemistry{<:Any, <:FlatSediment}}, Δt, χ)
    workgroup, worksize = work_layout(model.grid, :xyz)
    arch = model.architecture
    barrier = device_event(arch)
    step_field_kernel! = ab2_step_field!(device(arch), workgroup, worksize)
    model_fields = prognostic_fields(model)
    events = []

    for (i, field) in enumerate(model_fields)

        field_event = step_field_kernel!(field, Δt, χ,
                                         model.timestepper.Gⁿ[i],
                                         model.timestepper.G⁻[i],
                                         dependencies = device_event(arch))

        push!(events, field_event)

        # TODO: function tracer_index(model, field_index) = field_index - 3, etc...
        tracer_index = Val(i - 3) # assumption

        implicit_step!(field,
                       model.timestepper.implicit_solver,
                       model.closure,
                       model.diffusivity_fields,
                       tracer_index,
                       model.clock,
                       Δt,
                       dependencies = field_event)
    end

    sediment = model.biogeochemistry.sediment_model

    for (i, field) in enumerate(sediment_fields(model.biogeochemistry.sediment_model))
        field_event = launch!(arch, grid, :xy, ab2_step_sediment_field!, field, Δt, χ, sediment.tendencies.Gⁿ[i], sediment.tendencies.G⁻[i], dependencies=device_event(arch))

        push!(events, field_event)
    end

    wait(device(model.architecture), MultiEvent(Tuple(events)))

    return nothing
end

@kernel function ab2_step_sediment_field!(u, Δt, χ, Gⁿ, G⁻)
    i, j = @index(Global, NTuple)

    T = eltype(u)
    one_point_five = convert(T, 1.5)
    oh_point_five = convert(T, 0.5)

    @inbounds u[i, j, 1] += Δt * ((one_point_five + χ) * Gⁿ[i, j, 1] - (oh_point_five + χ) * G⁻[i, j, 1])
end
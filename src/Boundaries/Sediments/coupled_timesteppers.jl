using Oceananigans: NonhydrostaticModel, prognostic_fields
using OceanBioME: ContinuousFormBiogeochemistry
using OceanBioME.Boundaries.Sediments: AbstractSediment
using Oceananigans.TimeSteppers: ab2_step_field!, rk3_substep_field!, stage_Δt
using Oceananigans.Utils: work_layout, launch!
using Oceananigans.TurbulenceClosures: implicit_step!

import Oceananigans.TimeSteppers: ab2_step!, rk3_substep!

@inline function ab2_step!(model::NonhydrostaticModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:ContinuousFormBiogeochemistry{<:Any, <:FlatSediment}}, Δt, χ)
    workgroup, worksize = work_layout(model.grid, :xyz)
    arch = model.architecture
    step_field_kernel! = ab2_step_field!(device(arch), workgroup, worksize)
    model_fields = prognostic_fields(model)

    for (i, field) in enumerate(model_fields)

        step_field_kernel!(field, Δt, χ,
                           model.timestepper.Gⁿ[i],
                           model.timestepper.G⁻[i],
                           dependencies = barrier)

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

    for (i, field) in enumerate(sediment_fields(sediment))
        launch!(arch, model.grid, :xy, ab2_step_flat_field!, 
                field, Δt, χ, 
                sediment.tendencies.Gⁿ[i], 
                sediment.tendencies.G⁻[i])

    end

    return nothing
end

@kernel function ab2_step_flat_field!(u, Δt, χ, Gⁿ, G⁻)
    i, j = @index(Global, NTuple)

    T = eltype(u)
    one_point_five = convert(T, 1.5)
    oh_point_five = convert(T, 0.5)

    @inbounds u[i, j, 1] += Δt * ((one_point_five + χ) * Gⁿ[i, j, 1] - (oh_point_five + χ) * G⁻[i, j, 1])
end


function rk3_substep!(model::NonhydrostaticModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:ContinuousFormBiogeochemistry{<:Any, <:FlatSediment}}, Δt, γⁿ, ζⁿ)
    workgroup, worksize = work_layout(model.grid, :xyz)
    arch = model.architecture
    substep_field_kernel! = rk3_substep_field!(device(arch), workgroup, worksize)
    model_fields = prognostic_fields(model)

    for (i, field) in enumerate(model_fields)
        substep_field_kernel!(field, Δt, γⁿ, ζⁿ,
                              model.timestepper.Gⁿ[i],
                              model.timestepper.G⁻[i])

        # TODO: function tracer_index(model, field_index) = field_index - 3, etc...
        tracer_index = Val(i - 3) # assumption

        implicit_step!(field,
                       model.timestepper.implicit_solver,
                       model.closure,
                       model.diffusivity_fields,
                       tracer_index,
                       model.clock,
                       stage_Δt(Δt, γⁿ, ζⁿ))
    end

    sediment = model.biogeochemistry.sediment_model

    for (i, field) in enumerate(sediment_fields(sediment))
        launch!(arch, model.grid, :xy, rk3_step_flat_field!, 
                field, Δt, γⁿ, ζⁿ,
                sediment.tendencies.Gⁿ[i], 
                sediment.tendencies.G⁻[i])
    end


    return nothing
end

@kernel function rk3_step_flat_field!(U, Δt, γⁿ, ζⁿ, Gⁿ, G⁻)
    i, j = @index(Global, NTuple)

    @inbounds U[i, j, 1] += Δt * (γⁿ * Gⁿ[i, j, 1] + ζⁿ * G⁻[i, j, 1])
end

@kernel function rk3_step_flat_field!(U, Δt, γ¹, ::Nothing, G¹, G⁰)
    i, j = @index(Global, NTuple)

    @inbounds U[i, j, 1] += Δt * γ¹ * G¹[i, j, 1]
end
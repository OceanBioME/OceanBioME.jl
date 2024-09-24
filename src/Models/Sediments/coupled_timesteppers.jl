using Oceananigans: NonhydrostaticModel, prognostic_fields, HydrostaticFreeSurfaceModel
using OceanBioME.Models.Sediments: AbstractSediment
using Oceananigans.TimeSteppers: ab2_step_field!, rk3_substep_field!, stage_Δt
using Oceananigans.Utils: work_layout, launch!
using Oceananigans.TurbulenceClosures: implicit_step!
using Oceananigans.Models.HydrostaticFreeSurfaceModels: local_ab2_step!, ab2_step_free_surface!
using Oceananigans.Architectures: AbstractArchitecture

import Oceananigans.TimeSteppers: ab2_step!, rk3_substep!

const BGC_WITH_FLAT_SEDIMENT = Union{<:DiscreteBiogeochemistry{<:Any, <:Any, <:FlatSediment},
                                     <:ContinuousBiogeochemistry{<:Any, <:Any, <:FlatSediment}}

# This is definitly type piracy
@inline function ab2_step!(model::NonhydrostaticModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, BGC_WITH_FLAT_SEDIMENT}, Δt)
    workgroup, worksize = work_layout(model.grid, :xyz)
    arch = model.architecture
    step_field_kernel! = ab2_step_field!(device(arch), workgroup, worksize)
    model_fields = prognostic_fields(model)
    χ = model.timestepper.χ

    for (i, field) in enumerate(model_fields)

        step_field_kernel!(field, Δt, χ,
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
                       Δt)
    end

    sediment = model.biogeochemistry.sediment

    for (i, field) in enumerate(sediment_fields(sediment))
        launch!(arch, model.grid, :xy, ab2_step_flat_field!, 
                field, Δt, χ, 
                sediment.tendencies.Gⁿ[i], 
                sediment.tendencies.G⁻[i])

    end

    return nothing
end

@inline function ab2_step!(model::HydrostaticFreeSurfaceModel{<:Any, <:Any, <:AbstractArchitecture, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, BGC_WITH_FLAT_SEDIMENT}, Δt)
    χ = model.timestepper.χ

    # Step locally velocity and tracers
    @apply_regionally local_ab2_step!(model, Δt, χ)

    # blocking step for implicit free surface, non blocking for explicit
    ab2_step_free_surface!(model.free_surface, model, Δt, χ)

    sediment = model.biogeochemistry.sediment
    arch = model.architecture
    

    for (i, field) in enumerate(sediment_fields(sediment))
        launch!(arch, model.grid, :xy, ab2_step_flat_field!, 
                field, Δt, χ, 
                sediment.tendencies.Gⁿ[i], 
                sediment.tendencies.G⁻[i])

    end
end

@kernel function ab2_step_flat_field!(u, Δt, χ, Gⁿ, G⁻)
    i, j = @index(Global, NTuple)

    T = eltype(u)
    one_point_five = convert(T, 1.5)
    oh_point_five = convert(T, 0.5)

    @inbounds u[i, j, 1] += Δt * ((one_point_five + χ) * Gⁿ[i, j, 1] - (oh_point_five + χ) * G⁻[i, j, 1])
end

function rk3_substep!(model::NonhydrostaticModel{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, BGC_WITH_FLAT_SEDIMENT}, Δt, γⁿ, ζⁿ)
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

    sediment = model.biogeochemistry.sediment

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
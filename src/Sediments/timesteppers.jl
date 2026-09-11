using Oceananigans.Architectures: architecture
using Oceananigans.TimeSteppers: QuasiAdamsBashforth2TimeStepper, RungeKutta3TimeStepper,
                                 time_step!, tick!
using Oceananigans.Utils: launch!

import Oceananigans.TimeSteppers: ab2_step!, rk3_substep!,
                                  cache_previous_tendencies!, compute_flux_bc_tendencies!

const SUPPORTED_TIMESTEPPERS = Union{QuasiAdamsBashforth2TimeStepper, RungeKutta3TimeStepper}

function step_sediment!(sediment, model, parent_timestepper::SUPPORTED_TIMESTEPPERS)
    substep_sediment!(sediment, model, parent_timestepper)

    cache_previous_tendencies!(sediment)

    tick!(sediment.clock, model.clock.last_stage_Δt)

    return nothing
end

function step_sediment!(sediment, model, parent_timestepper)
    @warn "$(nameof(typeof(parent_timestepper))) is not supported for sediment models, " *
          "so the sediment is being stepped independently of $(nameof(typeof(model))), " *
          "which is neither mass conserving nor higher than first order accurate" maxlog = 1

    time_step!(sediment, model.clock.last_stage_Δt)

    return nothing
end

function substep_sediment!(sediment, model, parent_timestepper::QuasiAdamsBashforth2TimeStepper)
    ab2_step!(sediment, model.clock.last_Δt, parent_timestepper.χ, nothing)

    return nothing
end

function substep_sediment!(sediment, model, parent_timestepper::RungeKutta3TimeStepper)
    γ¹ = parent_timestepper.γ¹
    γ² = parent_timestepper.γ²
    γ³ = parent_timestepper.γ³
    ζ² = parent_timestepper.ζ²
    ζ³ = parent_timestepper.ζ³

    stage_Δt = model.clock.last_stage_Δt

    if model.clock.stage == 2
        rk3_substep!(sediment, stage_Δt / γ¹, γ¹, nothing, nothing)
    elseif model.clock.stage == 3
        rk3_substep!(sediment, stage_Δt / (γ² + ζ²), γ², ζ², nothing)
    else
        rk3_substep!(sediment, model.clock.last_Δt, γ³, ζ³, nothing)
    end

    return nothing
end

# AB2 methods

ab2_step!(model::BiogeochemicalSediment, Δt, callbacks) =
    ab2_step!(model, Δt, model.timestepper.χ, callbacks)

function ab2_step!(model::BiogeochemicalSediment, Δt, χ, callbacks)
    grid = model.grid
    arch = architecture(grid)
    model_fields = prognostic_fields(model)

    for (i, field) in enumerate(model_fields)
        kernel_args = (field, Δt, χ, model.timestepper.Gⁿ[i], model.timestepper.G⁻[i])
        launch!(arch, grid, :xy, ab2_step_flat_field!, kernel_args...; exclude_periphery=true)
    end

    return nothing
end

@kernel function ab2_step_flat_field!(u, Δt, χ, Gⁿ, G⁻)
    i, j = @index(Global, NTuple)

    FT = typeof(χ)
    Δt = convert(FT, Δt)
    one_point_five = convert(FT, 1.5)
    oh_point_five  = convert(FT, 0.5)
    not_euler = χ != convert(FT, -0.5) # use to prevent corruption by leftover NaNs in G⁻

    @inbounds begin
        Gu = (one_point_five + χ) * Gⁿ[i, j] - (oh_point_five + χ) * G⁻[i, j] * not_euler
        u[i, j, 1] += Δt * Gu
    end
end

# RK3 methods

function rk3_substep!(model::BiogeochemicalSediment, Δt, γⁿ, ζⁿ, callbacks)
    grid = model.grid
    arch = architecture(grid)
    model_fields = prognostic_fields(model)

    for (i, field) in enumerate(model_fields)
        kernel_args = (field, Δt, γⁿ, ζⁿ, model.timestepper.Gⁿ[i], model.timestepper.G⁻[i])
        launch!(arch, grid, :xy, rk3_substep_flat_field!, kernel_args...; exclude_periphery=true)
    end

    return nothing
end

@kernel function rk3_substep_flat_field!(U, Δt, γⁿ::FT, ζⁿ, Gⁿ, G⁻) where FT
    i, j = @index(Global, NTuple)

    @inbounds begin
        U[i, j, 1] += convert(FT, Δt) * (γⁿ * Gⁿ[i, j] + ζⁿ * G⁻[i, j])
    end
end

@kernel function rk3_substep_flat_field!(U, Δt, γ¹::FT, ::Nothing, G¹, G⁰) where FT
    i, j = @index(Global, NTuple)

    @inbounds begin
        U[i, j, 1] += convert(FT, Δt) * γ¹ * G¹[i, j]
    end
end

# store tendencies

""" Store source terms for `u`, `v`, and `w`. """
@kernel function store_flat_field_tendencies!(G⁻, G⁰)
    i, j = @index(Global, NTuple)
    @inbounds G⁻[i, j, 1] = G⁰[i, j, 1]
end

""" Store previous source terms before updating them. """
function cache_previous_tendencies!(model::BiogeochemicalSediment)
    model_fields = prognostic_fields(model)

    for field_name in keys(model_fields)
        launch!(architecture(model.grid), model.grid, :xy, store_flat_field_tendencies!,
                model.timestepper.G⁻[field_name],
                model.timestepper.Gⁿ[field_name])
    end

    return nothing
end

compute_flux_bc_tendencies!(model::BiogeochemicalSediment) = nothing

deprecate_sediment_timestepper(::Nothing) = nothing

deprecate_sediment_timestepper(timestepper) =
    @warn "The `timestepper` keyword argument to sediment models is deprecated and ignored: " *
          "the sediment is stepped with the time stepping coefficients of the model it is " *
          "coupled to, so `timestepper = :$(timestepper)` has no effect" maxlog = 1

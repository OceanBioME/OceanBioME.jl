using Oceananigans.Architectures: architecture
using Oceananigans.TimeSteppers: QuasiAdamsBashforth2TimeStepper, RungeKutta3TimeStepper
using Oceananigans.Utils: launch!

import Oceananigans.TimeSteppers: ab2_step!, rk3_substep!,
                                  store_tendencies!

const VALID_TIMESTEPPERS = Union{<:QuasiAdamsBashforth2TimeStepper, <:RungeKutta3TimeStepper}

validate_sediment_timestepper(timestepper) = throw(ArgumentError("$(typeof(timestepper)) is not configured for sediment models"))
validate_sediment_timestepper(::VALID_TIMESTEPPERS) = nothing 

# AB2 methods

function ab2_step!(model::BiogeochemicalSediment, Δt)
    grid = model.grid
    arch = architecture(grid)
    model_fields = prognostic_fields(model)
    χ = model.timestepper.χ

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
        u[i, j] += Δt * Gu
    end
end

# RK3 methods

function rk3_substep!(model::BiogeochemicalSediment, Δt, γⁿ, ζⁿ)
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
        U[i, j] += convert(FT, Δt) * (γⁿ * Gⁿ[i, j] + ζⁿ * G⁻[i, j])
    end
end

@kernel function rk3_substep_flat_field!(U, Δt, γ¹::FT, ::Nothing, G¹, G⁰) where FT
    i, j = @index(Global, NTuple)

    @inbounds begin
        U[i, j] += convert(FT, Δt) * γ¹ * G¹[i, j]
    end
end

# store tendencies

""" Store source terms for `u`, `v`, and `w`. """
@kernel function store_flat_field_tendencies!(G⁻, G⁰)
    i, j = @index(Global, NTuple)
    @inbounds G⁻[i, j] = G⁰[i, j]
end

""" Store previous source terms before updating them. """
function store_tendencies!(model::BiogeochemicalSediment)
    model_fields = prognostic_fields(model)

    for field_name in keys(model_fields)
        launch!(architecture(model.grid), model.grid, :xy, store_flat_field_tendencies!,
                model.timestepper.G⁻[field_name],
                model.timestepper.Gⁿ[field_name])
    end

    return nothing
end
#=struct RungeKutta3TimeStepper{FT, TG}
    γ¹ :: FT
    γ² :: FT
    γ³ :: FT
    ζ² :: FT
    ζ³ :: FT
    Gⁿ :: TG
    G⁻ :: TG
end

function RungeKutta3TimeStepper(variables)
    γ¹ = 8 // 15
    γ² = 5 // 12
    γ³ = 3 // 4

    ζ² = -17 // 60
    ζ³ = -5 // 12

    Gⁿ = NamedTuple{variables}([CenterField(BoxModelGrid) for var in eachindex(variables)])
    G⁻ = NamedTuple{variables}([CenterField(BoxModelGrid) for var in eachindex(variables)])

    FT = eltype(Gⁿ[1])
    TG = typeof(Gⁿ)

    return RungeKutta3TimeStepper{FT, TG}(γ¹, γ², γ³, ζ², ζ³, Gⁿ, G⁻)
end

function time_step!(model::BoxModel{<:Any, <:Any, <:Any, <:RungeKutta3TimeStepper, <:Any}, Δt)
    Δt == 0 && error("Δt can not be zero")

    γ¹ = model.timestepper.γ¹
    γ² = model.timestepper.γ²
    γ³ = model.timestepper.γ³

    ζ² = model.timestepper.ζ²
    ζ³ = model.timestepper.ζ³

    first_stage_Δt  = γ¹ * Δt
    second_stage_Δt = (γ² + ζ²) * Δt
    third_stage_Δt  = (γ³ + ζ³) * Δt

    # first stage
    calculate_tendencies!(model)
    rk3_substep!(model, Δt, γ¹, nothing)

    tick!(model.clock, first_stage_Δt; stage=true)
    store_tendencies!(model)

    update_boxmodel_state!(model)

    # second stage
    calculate_tendencies!(model)
    rk3_substep!(model, Δt, γ², ζ²)

    tick!(model.clock, second_stage_Δt; stage=true)
    store_tendencies!(model)

    update_boxmodel_state!(model)

    # third stage
    calculate_tendencies!(model)
    rk3_substep!(model, Δt, γ³, ζ³)

    tick!(model.clock, third_stage_Δt)
    store_tendencies!(model)

    update_boxmodel_state!(model)

    return nothing
end

function rk3_substep!(model, Δt, γⁿ, ζⁿ)
    @inbounds for (i, field) in enumerate(model.fields)
        field[1] += Δt * (γⁿ * model.timestepper.Gⁿ[i][1] + ζⁿ * model.timestepper.G⁻[i][1])
    end
end

function rk3_substep!(model, Δt, γⁿ, ::Nothing)
    for (i, field) in enumerate(model.fields)
        @inbounds field[1] += Δt * γⁿ * model.timestepper.Gⁿ[i][1]
    end
end

function store_tendencies!(model)
    @inbounds for (i, Gⁿ) in enumerate(model.timestepper.Gⁿ)
        model.timestepper.G⁻[i][1] = Gⁿ[i]
    end
end

summary(::RungeKutta3TimeStepper{FT, TG}) where {FT, TG} = string("Runge-Kutta 3 Timetepper")
show(io::IO, model::RungeKutta3TimeStepper{FT, TG}) where {FT, TG} = print(io, summary(model))
=#
using Oceananigans.Architectures: device
using Oceananigans.Biogeochemistry: update_tendencies!
using Oceananigans.TimeSteppers: rk3_substep_field!, store_field_tendencies!, RungeKutta3TimeStepper, QuasiAdamsBashforth2TimeStepper
using Oceananigans.Utils: work_layout, launch!

import Oceananigans.TimeSteppers: rk3_substep!, store_tendencies!, compute_tendencies!

function BoxModelTimeStepper(name::Symbol, grid, tracers)
    fullname = Symbol(name, :TimeStepper)

    Gⁿ = TracerFields(tracers, grid)
    G⁻ = TracerFields(tracers, grid)
    
    return @eval $fullname($grid, $tracers; Gⁿ = $Gⁿ, G⁻ = $G⁻)
end

""" Store previous source terms before updating them. """
function store_tendencies!(model::BoxModel)
    model_fields = prognostic_fields(model)

    for field_name in keys(model_fields)
        launch!(architecture(model), model.grid, :xyz, store_field_tendencies!,
                model.timestepper.G⁻[field_name],
                model.grid,
                model.timestepper.Gⁿ[field_name])
    end

    return nothing
end

function compute_tendencies!(model::BoxModel, callbacks)
    Gⁿ = model.timestepper.Gⁿ

    model_fields = @inbounds [field[1] for field in model.fields]

    @inbounds for tracer in required_biogeochemical_tracers(model.biogeochemistry)
        if !(tracer == :T)
            getproperty(Gⁿ, tracer)[1] = model.biogeochemistry(Val(tracer), 0.0, 0.0, 0.0, model.clock.time, model_fields...) + model.forcing[tracer](model.clock.time, model_fields...)
        end
    end

    @inbounds for variable in required_biogeochemical_auxiliary_fields(model.biogeochemistry)
        if !(variable == :PAR)
            getproperty(Gⁿ, variable)[1] = model.forcing[variable](model.clock.time, model_fields...)
        end
    end

    for callback in callbacks
        callback.callsite isa TendencyCallsite && callback(model)
    end

    update_tendencies!(model.biogeochemistry, model)

    return nothing
end

function rk3_substep!(model::BoxModel, Δt, γⁿ, ζⁿ)

    workgroup, worksize = work_layout(model.grid, :xyz)
    substep_field_kernel! = rk3_substep_field!(device(architecture(model)), workgroup, worksize)
    model_fields = prognostic_fields(model)

    for (i, field) in enumerate(model_fields)
        substep_field_kernel!(field, Δt, γⁿ, ζⁿ,
                              model.timestepper.Gⁿ[i],
                              model.timestepper.G⁻[i])
    end

    return nothing
end

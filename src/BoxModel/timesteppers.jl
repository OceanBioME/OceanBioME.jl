struct RungeKutta3TimeStepper{FT, TG}
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

    Gⁿ = StructArray([NamedTuple{variables}(zeros(length(variables)))])
    G⁻ = StructArray([NamedTuple{variables}(zeros(length(variables)))])

    FT = typeof(Gⁿ[1][1])
    TG = typeof(Gⁿ)

    return RungeKutta3TimeStepper{FT, TG}(γ¹, γ², γ³, ζ², ζ³, Gⁿ, G⁻)
end

function time_step!(model::BoxModel{<:Any, <:Any, <:Any, <:Any, <:RungeKutta3TimeStepper, <:Any}, Δt)
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
    @inbounds for tracer in keys(model.values[1])
        getproperty(model.values, tracer) .+= Δt * (γⁿ * getproperty(model.timestepper.Gⁿ, tracer) + ζⁿ * getproperty(model.timestepper.G⁻, tracer))
    end
end

function rk3_substep!(model, Δt, γⁿ, ::Nothing)
    @inbounds for tracer in keys(model.values[1])
        getproperty(model.values, tracer) .+= Δt * γⁿ * getproperty(model.timestepper.Gⁿ, tracer)
    end
end

function store_tendencies!(model)
    @inbounds for tracer in keys(model.values[1])
        getproperty(model.timestepper.G⁻, tracer) .= getproperty(model.timestepper.Gⁿ, tracer)
    end
end
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

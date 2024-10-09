using Oceananigans.Architectures: device
using Oceananigans.Biogeochemistry: update_tendencies!, biogeochemical_auxiliary_fields, AbstractBiogeochemistry, AbstractContinuousFormBiogeochemistry
using Oceananigans.Grids: nodes, Center
using Oceananigans.TimeSteppers: rk3_substep_field!, store_field_tendencies!, RungeKutta3TimeStepper, QuasiAdamsBashforth2TimeStepper
using Oceananigans.Utils: work_layout, launch!
using Oceananigans: TendencyCallsite

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

    @inbounds for field_name in keys(model_fields)
        model.timestepper.G⁻[field_name][1, 1, 1] = model.timestepper.Gⁿ[field_name][1, 1, 1]
    end

    return nothing
end

function compute_tendencies!(model::BoxModel, callbacks)
    Gⁿ = model.timestepper.Gⁿ

    @inbounds for (i, field) in enumerate(model.fields)
        model.field_values[i] = field[1, 1, 1]
    end

    @inbounds for (i, field) in enumerate(values(biogeochemical_auxiliary_fields(model.biogeochemistry)))
        model.field_values[i+length(model.fields)] = field[1, 1, 1]
    end

    for tracer in required_biogeochemical_tracers(model.biogeochemistry)
        forcing = @inbounds model.forcing[tracer]
        
        @inbounds Gⁿ[tracer][1, 1, 1] = tracer_tendency(Val(tracer), model.biogeochemistry, forcing, model.clock, model.fields, model.field_values, model.grid)
    end

    for callback in callbacks
        callback.callsite isa TendencyCallsite && callback(model)
    end

    update_tendencies!(model.biogeochemistry, model)

    return nothing
end

@inline boxmodel_xyz(nodes, grid) = map(n->boxmodel_coordinate(n, grid), nodes)
@inline boxmodel_coordinate(::Nothing, grid) = zero(grid)
@inline boxmodel_coordinate(nodes, grid) = @inbounds nodes[1]

@inline tracer_tendency(val_name, biogeochemistry::AbstractContinuousFormBiogeochemistry, forcing, clock, model_fields, model_field_values, grid) =
    biogeochemistry(val_name, boxmodel_xyz(nodes(grid, Center(), Center(), Center()), grid)..., clock.time, model_field_values...) + forcing(time, model_fields...)

@inline tracer_tendency(val_name, biogeochemistry::AbstractBiogeochemistry, forcing, clock, model_fields, model_field_values, grid) =
    biogeochemistry(1, 1, 1, grid, val_name, clock, model_fields) + forcing(clock, model_fields)

function rk3_substep!(model::BoxModel, Δt, γⁿ, ζⁿ)
    model_fields = prognostic_fields(model)

    for (i, field) in enumerate(model_fields)
        Gⁿ = @inbounds model.timestepper.Gⁿ[i]
        G⁻ = @inbounds model.timestepper.G⁻[i]

        rk3_substep!(field, Δt, γⁿ, ζⁿ, Gⁿ, G⁻)
    end

    return nothing
end

function rk3_substep!(U, Δt, γⁿ::FT, ζⁿ, Gⁿ, G⁻) where FT
    @inbounds begin
        U[1, 1, 1] += convert(FT, Δt) * (γⁿ * Gⁿ[1, 1, 1] + ζⁿ * G⁻[1, 1, 1])
    end
    
    return nothing
end

function rk3_substep!(U, Δt, γ¹::FT, ::Nothing, G¹, G⁰) where FT
    @inbounds begin
        U[1, 1, 1] += convert(FT, Δt) * γ¹ * G¹[1, 1, 1]
    end

    return nothing
end

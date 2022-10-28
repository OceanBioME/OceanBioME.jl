using Oceananigans.Fields: TendencyFields
using Oceananigans.TimeSteppers: rk3_substep_field!
using OceanBioME.Oceananigans.Fields: isacolumn

import Oceananigans.TimeSteppers: RungeKutta3TimeStepper, rk3_substep!, stage_Δt

function RungeKutta3TimeStepper(grid, tracers, forced_auxiliary_fields;
    implicit_solver::TI = nothing,
    Gⁿ::TG = merge(TendencyFields(grid, tracers), AuxiliaryFields(forced_auxiliary_fields, grid)),
    G⁻ = merge(TendencyFields(grid, tracers), AuxiliaryFields(forced_auxiliary_fields, grid))) where {TI, TG}

    !isnothing(implicit_solver) &&
    @warn("Implicit-explicit time-stepping with RungeKutta3TimeStepper is not tested. " * 
    "\n implicit_solver: $(typeof(implicit_solver))")

    γ¹ = 8 // 15
    γ² = 5 // 12
    γ³ = 3 // 4

    ζ² = -17 // 60
    ζ³ = -5 // 12

    FT = eltype(grid)

    return RungeKutta3TimeStepper{FT, TG, TI}(γ¹, γ², γ³, ζ², ζ³, Gⁿ, G⁻, implicit_solver)
end


function rk3_substep!(model, Δt, γⁿ, ζⁿ)

    workgroup, worksize = work_layout(model.grid, :xyz)
    arch = model.architecture
    barrier = Event(device(arch))
    substep_field_kernel! = rk3_substep_field!(device(arch), workgroup, worksize)
    substep_column_field_kernel! = rk3_substep_field!(device(arch), workgroup, (worksize[1], worksize[2], 1))
    model_fields = prognostic_fields(model)
    events = []

    for (i, field) in enumerate(model_fields)

        field_event = isacolumn(field) ? substep_column_field_kernel!(field, Δt, γⁿ, ζⁿ,
                                                                            model.timestepper.Gⁿ[i],
                                                                            model.timestepper.G⁻[i],
                                                                            dependencies=barrier) : substep_field_kernel!(field, Δt, γⁿ, ζⁿ,
                                                                                                                                                model.timestepper.Gⁿ[i],
                                                                                                                                                model.timestepper.G⁻[i],
                                                                                                                                                dependencies=barrier)

        # TODO: function tracer_index(model, field_index) = field_index - 3, etc...
        tracer_index = Val(i - 3) # assumption

        implicit_step!(field,
                       model.timestepper.implicit_solver,
                       model.closure,
                       model.diffusivity_fields,
                       tracer_index,
                       model.clock,
                       stage_Δt(Δt, γⁿ, ζⁿ),
                       dependencies = field_event)

        push!(events, field_event)
    end

    wait(device(arch), MultiEvent(Tuple(events)))

    return nothing
end

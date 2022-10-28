using Oceananigans.Fields: TendencyFields
using Oceananigans.TimeSteppers: ab2_step_field!

import Oceananigans.TimeSteppers: QuasiAdamsBashforth2TimeStepper, ab2_step!

function QuasiAdamsBashforth2TimeStepper(grid, tracers, forced_auxiliary_fields,
        χ = 0.1;
        implicit_solver::IT = nothing,
        Gⁿ = merge(TendencyFields(grid, tracers), AuxiliaryFields(forced_auxiliary_fields, grid)),
        G⁻ = merge(TendencyFields(grid, tracers), AuxiliaryFields(forced_auxiliary_fields, grid))) where IT

    FT = eltype(grid)
    GT = typeof(Gⁿ)

    return QuasiAdamsBashforth2TimeStepper{FT, GT, IT}(χ, Inf, Gⁿ, G⁻, implicit_solver)
end

""" Generic implementation. """
function ab2_step!(model, Δt, χ)

    workgroup, worksize = work_layout(model.grid, :xyz)
    arch = model.architecture
    barrier = device_event(arch)
    step_field_kernel! = ab2_step_field!(device(arch), workgroup, worksize)
    step_column_field_kernel! = ab2_step_field!(device(arch), workgroup, (worksize[1], worksize[2], 1))
    model_fields = prognostic_fields(model)
    events = []

    for (i, field) in enumerate(model_fields)
        # this is faster than ifelse
        field_event = isacolumn(field) ? step_column_field_kernel!(field, Δt, χ,
                                                                    model.timestepper.Gⁿ[i],
                                                                    model.timestepper.G⁻[i],
                                                                    dependencies = device_event(arch)) : step_field_kernel!(field, Δt, χ,
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

    wait(device(model.architecture), MultiEvent(Tuple(events)))

    return nothing
end
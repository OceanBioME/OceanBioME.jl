import Oceananigans.TimeSteppers: store_tendencies!, store_field_tendencies!

""" Store previous source terms before updating them. """
function store_tendencies!(model)

    workgroup, worksize = work_layout(model.grid, :xyz)
    arch = model.architecture
    store_kernel! = store_field_tendencies!(device(arch), workgroup, worksize)
    store_column_kernel! = store_field_tendencies!(device(arch), workgroup, (worksize[1], worksize[2], 1))

    model_fields = prognostic_fields(model)

    events = []

    for field_name in keys(model_fields)
        G⁻ = model.timestepper.G⁻[field_name]
        field_event = isacolumn(G⁻) ? store_column_kernel!(G⁻,
                                                                    model.grid,
                                                                    model.timestepper.Gⁿ[field_name],
                                                                    dependencies = device_event(arch)) : store_kernel!(model.timestepper.G⁻[field_name],
                                                                                                                                            model.grid,
                                                                                                                                            model.timestepper.Gⁿ[field_name],
                                                                                                                                            dependencies = device_event(arch))

        push!(events, field_event)
    end

    wait(device(model.architecture), MultiEvent(Tuple(events)))

    return nothing
end
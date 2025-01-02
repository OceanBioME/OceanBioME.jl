using Oceananigans.Biogeochemistry: extract_biogeochemical_fields
using Oceananigans.Fields: Center
using Oceananigans.Grids: xnode, ynode

function compute_sediment_tendencies!(model)
    field_names = required_sediment_fields(model)

    sediment_fields = model.fields
    tracked_fields = model.tracked_fields

    grid = model.grid
    clock = model.clock

    biogeochemistry = model.biogeochemistry

    Gⁿ = model.timestepper.Gⁿ

    arch = architecture(grid)

    for field_name in field_names
        G = Gⁿ[field_name]

        args = (Val(field_name), biogeochemistry, sediment_fields, tracked_fields, clock)

        launch!(arch, grid, :xy, compute_sediment_tendency!, G, grid, args)
    end

    return nothing
end

@kernel function compute_sediment_tendency!(G, grid, args)
    i, j = @index(Global, NTuple)

    @inbounds G[i, j] = sediment_tendencies(i, j, grid, args...)
end

@inline sediment_tendencies(i, j, grid, val_name, biogeochemistry, fields, tracked_fields, clock) =
    biogeochemistry(i, j, grid, val_name, clock, fields, tracked_fields)

@inline function sediment_tendencies(i, j, grid, val_name, biogeochemistry::ACFSBGC, fields, tracked_fields, clock)
    model_fields = merge(fields, tracked_fields)

    field_values = extract_biogeochemical_fields(i, j, 1, grid, model_fields, keys(model_fields))

    x = xnode(i, j, 1, grid, Center(), Center(), Center())
    y = ynode(i, j, 1, grid, Center(), Center(), Center())

    return biogeochemistry(val_name, x, y, clock.time, field_values...)
end

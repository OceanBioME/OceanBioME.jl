using Oceananigans.Operators: Δzᶜᶜᶠ

function update_tendencies!(bgc, sediment_model::BiogeochemicalSediment, model)
    coupled_fields = coupled_tracers(sediment_model)
    sediment_fields = fields(sediment_model)
    tracked_fields = sediment_model.tracked_fields

    grid = model.grid
    clock = model.clock

    biogeochemistry = sediment_model.biogeochemistry

    Gⁿ = model.timestepper.Gⁿ

    bottom_indices = sediment_model.bottom_indices

    arch = architecture(grid)

    for field_name in coupled_fields
        G = Gⁿ[field_name]

        args = (Val(field_name), biogeochemistry, sediment_fields, tracked_fields, clock)

        launch!(arch, grid, :xy, update_coupled_tendency!, G, grid, bottom_indices, args)
    end

    return nothing
end

@kernel function update_coupled_tendency!(G, grid, bottom_indices, args)
    i, j = @index(Global, NTuple)

    @inbounds begin
        k = bottom_indices[i, j]

        @inbounds G[i, j, k] += sediment_tendencies(i, j, grid, args...) / Δzᶜᶜᶠ(i, j, k, grid)
    end
end
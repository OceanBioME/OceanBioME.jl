using Oceananigans.Advection: advective_tracer_flux_z, FluxFormAdvection
using Oceananigans.Architectures: architecture
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Models: total_velocities, AbstractModel
using Oceananigans.Operators: Azᶜᶜᶠ
using Oceananigans.Utils: launch!

function update_tracked_fields!(sediment, model)
    grid = model.grid
    arch = architecture(grid)

    model_fields = prognostic_fields(model)

    bottom_indices = sediment.bottom_indices

    # tracked tracers
    field_names = required_tracers(sediment)

    for field_name in field_names
        source = model_fields[field_name]
        destination = sediment.tracked_fields[field_name]

        launch!(arch, grid, :xy, copy_to_sediment!, source, destination, bottom_indices)
    end

    # tracked fluxs
    field_names = sinking_fluxs(sediment)

    for field_name in field_names
        source = model_fields[field_name]
        advection = vertical_advection_scheme(model, field_name)
        w = biogeochemical_drift_velocity(model.biogeochemistry, Val(field_name)).w
        destination = sediment.tracked_fields[field_name]

        launch!(arch, grid, :xy, compute_sinking_flux!, destination, source, advection, w, bottom_indices, grid)
    end

    return nothing
end

# tracer fields

@kernel function copy_to_sediment!(source, destination, bottom_indices)
    i, j = @index(Global, NTuple)

    @inbounds begin
        k = bottom_indices[i, j, 1]

        destination[i, j, 1] = source[i, j, k]
    end
end

# fluxs

@inline vertical_advection_scheme(advection, name) = advection
@inline vertical_advection_scheme(advection::FluxFormAdvection, name) = advection.z
@inline vertical_advection_scheme(advection::NamedTuple, name) = advection[name]
@inline vertical_advection_scheme(model::AbstractModel, name) = vertical_advection_scheme(model.advection, name)

@inline sinking_flux(i, j, k, grid, advection, w, C) =
    - advective_tracer_flux_z(i, j, k, grid, advection, w, C) / Azᶜᶜᶠ(i, j, k, grid)

@kernel function compute_sinking_flux!(destination, source, advection, w, bottom_indices, grid)
    i, j = @index(Global, NTuple)

    @inbounds begin
        k = bottom_indices[i, j, 1]

        destination[i, j, 1] = sinking_flux(i, j, k, grid, advection, w, source)
    end
end
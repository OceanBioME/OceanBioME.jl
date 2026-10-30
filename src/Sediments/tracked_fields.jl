using Oceananigans.Advection: _advective_tracer_flux_z, FluxFormAdvection
using Oceananigans.Architectures: architecture
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Models: total_velocities, AbstractModel
using Oceananigans.Operators: Azᶜᶜᶠ
using Oceananigans.Utils: launch!

_search_floor_flux(::Nothing, name) = nothing

function _search_floor_flux(m, name)
    hasproperty(m, :floor_flux) || return nothing
    ff = getproperty(m, :floor_flux)
    hasproperty(ff, name) || return nothing
    return getproperty(ff, name)
end

function _search_floor_flux(mods::Tuple, name)
    for m in mods
        ff = _search_floor_flux(m, name)
        ff !== nothing && return ff
    end
    return nothing
end

function _implicit_floor_flux(model, name)
    bgc = model.biogeochemistry
    if hasproperty(bgc, :modifiers)
        ff = _search_floor_flux(bgc.modifiers, name)
        ff !== nothing && return ff
    end
    if hasproperty(bgc, :underlying_biogeochemistry)
        det = bgc.underlying_biogeochemistry
        if hasproperty(det, :detritus) && hasproperty(det.detritus, :sinking)
            ff = _search_floor_flux(det.detritus.sinking, name)
            ff !== nothing && return ff
        end
    end
    return nothing
end

function update_tracked_fields!(sediment, model)
    grid = model.grid
    arch = architecture(grid)
    model_fields = fields(model)

    bottom_indices = sediment.bottom_indices

    # tracked tracers
    field_names = required_tracers(sediment)

    for field_name in field_names
        source = model_fields[field_name]
        destination = sediment.tracked_fields[field_name]

        launch!(arch, grid, :xy, copy_to_sediment!, source, destination, bottom_indices)
    end

    # tracked fluxes
    field_names = sinking_fluxes(sediment)

    for field_name in field_names
        destination = sediment.tracked_fields[field_name]

        ff = _implicit_floor_flux(model, field_name)

        if ff !== nothing
            launch!(arch, grid, :xy, _copy_floor_flux!, destination, ff)
        else
            source = model_fields[field_name]
            advection = vertical_advection_scheme(model, field_name)
            w = biogeochemical_drift_velocity(model.biogeochemistry, Val(field_name)).w
            launch!(arch, grid, :xy, compute_sinking_flux!, destination, source, advection, w, bottom_indices, grid)
        end
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

# fluxes

@inline vertical_advection_scheme(advection, name) = advection
@inline vertical_advection_scheme(advection::FluxFormAdvection, name) = advection.z
@inline vertical_advection_scheme(advection::NamedTuple, name) = advection[name]
@inline vertical_advection_scheme(model::AbstractModel, name) = vertical_advection_scheme(model.advection, name)

@inline sinking_flux(i, j, k, grid, advection, w, C) =
    - _advective_tracer_flux_z(i, j, k, grid, advection, w, C) / Azᶜᶜᶠ(i, j, k, grid)

@kernel function compute_sinking_flux!(destination, source, advection, w, bottom_indices, grid)
    i, j = @index(Global, NTuple)

    @inbounds begin
        k = bottom_indices[i, j, 1]

        destination[i, j, 1] = sinking_flux(i, j, k, grid, advection, w, source)
    end
end

@kernel function _copy_floor_flux!(destination, source)
    i, j = @index(Global, NTuple)
    @inbounds destination[i, j, 1] = source[i, j, 1]
end

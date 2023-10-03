using Oceananigans.Fields: ZFaceField, AbstractField
using Oceananigans.Forcings: maybe_constant_field
using Oceananigans.Grids: AbstractGrid

import Adapt: adapt_structure, adapt

function setup_velocity_fields(drift_speeds, grid::AbstractGrid, open_bottom; smoothing_distance = 2)
    drift_velocities = []
    for w in values(drift_speeds)
        if isa(values(w), Number)
            w_field = ZFaceField(grid)
            for k=1:grid.Nz 
                @inbounds w_field[:, :, k] .= - w * ifelse(open_bottom, 1.0, (1 - exp((1-k) / smoothing_distance)))
            end
            w = w_field
        elseif !isa(values(w), Tuple{AbstractField, AbstractField, AbstractField})
            error("Provided sinking speeds are an unsuitable value, must be a `NamedTuple` of scalar values (positive = sinking) with keys of tracer names, or `NamedTuple` of `VelocityFields`")
        end

        push!(drift_velocities, w)
    end

    return NamedTuple{keys(drift_speeds)}(drift_velocities)
end

setup_velocity_fields(drift_speeds, grid::BoxModelGrid, open_bottom) = drift_speeds

adapt_structure(to, velocities::NamedTuple{(:u, :v, :w), Tuple{AbstractField, AbstractField, AbstractField}}) = NamedTuple{(:u, :v, :w)}(adapt.(to, values(velocities)))

function show_sinking_velocities(sinking_velocities::NamedTuple{T, V}) where {T, V} 
    str = ""
    if length(T) == 0
        return str
    elseif length(T) == 1
        return "    └── $(T[1]): $(maximum_sinking(sinking_velocities[1])) to $(minimum_sinking(sinking_velocities[1])) m/s"
    else
        for idx in 1:length(T) - 1
            str *= "    ├── $(T[idx]): $(maximum_sinking(sinking_velocities[idx])) to $(minimum_sinking(sinking_velocities[idx])) m/s \n"
        end
        str *= "    └── $(T[end]): $(maximum_sinking(sinking_velocities[end])) to $(minimum_sinking(sinking_velocities[end])) m/s"

        return str
    end
end

maximum_sinking(velocity) = maximum(velocity)
minimum_sinking(velocity) = minimum(velocity)

maximum_sinking(velocity::NamedTuple) = maximum(velocity.w)
minimum_sinking(velocity::NamedTuple) = minimum(velocity.w)

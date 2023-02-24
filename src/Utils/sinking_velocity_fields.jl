using Oceananigans.Fields: ZFaceField, AbstractField
using Oceananigans.Forcings: maybe_constant_field
using Oceananigans.Grids: AbstractGrid

import Adapt: adapt_structure

function setup_velocity_fields(drift_speeds, grid::AbstractGrid, open_bottom; smoothing_distance = 2)
    drift_velocities = []
    for w in values(drift_speeds)
        u, v = maybe_constant_field.((0, 0))
        if isa(values(w), Number)
            w_field = ZFaceField(grid)
            for k=1:grid.Nz 
                @inbounds w_field[:, :, k] .= - w * ifelse(open_bottom, 1.0, (1 - exp((1-k) / smoothing_distance)))
            end
            w = w_field
        elseif !isa(values(w), Tuple{AbstractField, AbstractField, AbstractField})
            error("Provided sinking speeds are an unsuitable value, must be a `NamedTuple` of scalar values (positive = sinking) with keys of tracer names, or `NamedTuple` of `VelocityFields`")
        end

        push!(drift_velocities, (; u, v, w))
    end

    return NamedTuple{keys(drift_speeds)}(drift_velocities)
end

setup_velocity_fields(drift_speeds, grid::BoxModelGrid, open_bottom) = drift_speeds

adapt_structure(to, velocities::NamedTuple{(:u, :v, :w), Tuple{AbstractField, AbstractField, AbstractField}}) = NamedTuple{(:u, :v, :w)}(adapt_structure.(to, values(velocities)))

function show_sinking_velocities(sinking_velocities::NamedTuple{T, V}) where {T, V} 
    str = ""
    if length(T) == 1
        str = "    └── $(T[1]): $(-maximum(sinking_velocities[1].w)) to $(-minimum(sinking_velocities[1].w)) m/s"
    else
        for idx in 1:length(T) - 1
            str *= "    ├── $(T[idx]): $(-maximum(sinking_velocities[idx].w)) to $(-minimum(sinking_velocities[idx].w)) m/s \n"
        end
        str *= "    └── $(T[end]): $(-maximum(sinking_velocities[end].w)) to $(-minimum(sinking_velocities[end].w)) m/s"
    end
    return str
end
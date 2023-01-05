using Oceananigans.Fields: ZFaceField, AbstractField
using Oceananigans.Forcings: maybe_constant_field
using Oceananigans.Grids: AbstractGrid

function setup_velocity_fields(drift_speeds, grid::AbstractGrid, open_bottom)
    drift_velocities = []
    for w in values(drift_speeds)
        u, v = maybe_constant_field.((0, 0))
        if isa(values(w), Number)
            w_field = ZFaceField(grid)
            for k=1:grid.Nz 
                @inbounds w_field[:, :, k] .= - w * ifelse(open_bottom, 1.0, (1 - exp((1-k)/2)))
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
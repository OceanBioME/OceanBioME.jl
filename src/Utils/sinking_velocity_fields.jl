using Oceananigans.Fields: ZFaceField
using Oceananigans.Forcings: maybe_constant_field
using Oceananigans.Grids: AbstractGrid

function setup_velocity_fields(drift_speeds, grid::AbstractGrid, open_bottom)
    drift_velocities = []
    for w in values(drift_speeds)
        u, v = maybe_constant_field.((0, 0))
        if isa(values(w), Number)
            #if open_bottom
            #    w = maybe_constant_field(-w)
            #else
                w_field = ZFaceField(grid)
                for k=1:grid.Nz 
                    @inbounds w_field[:, :, k] .= - w * ifelse(open_bottom, 1.0, (1 - exp((1-k)/2)))
                end
                w = w_field
            #end
        end # otherwise we give up and rely on the user

        push!(drift_velocities, (; u, v, w))
    end

    return NamedTuple{keys(drift_speeds)}(drift_velocities)
end

setup_velocity_fields(drift_speeds, grid::BoxModelGrid, open_bottom) = drift_speeds
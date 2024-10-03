using Oceananigans.Operators: volume
using Oceananigans.Fields: fractional_indices, _interpolate, interpolator, Center
using Oceananigans.Grids: AbstractGrid, Flat, Bounded, Periodic

@inline get_node(::Bounded,  i, N) = min(max(i, 1), N)
@inline get_node(::Periodic, i, N) = ifelse(i < 1, N, ifelse(i > N, 1, i))

struct NearestPoint end

@inline function extract_tracer_values(::NearestPoint, particles, grid, fields, n)
    x = @inbounds particles.x[n]
    y = @inbounds particles.y[n]
    z = @inbounds particles.z[n]

    i, j, k = nearest_node(x, y, z, grid)

    field_names = required_tracers(particles)
    nf = length(field_names)

    field_values = ntuple(Val(nf)) do n
        fields[field_names[n]][i, j, k]
    end

    return field_values
end

@inline function nearest_node(x, y, z, grid::AbstractGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    # messy
    ii, jj, kk = fractional_indices((x, y, z), grid, ifelse(isa(TX(), Flat), nothing, Center()),
                                                     ifelse(isa(TY(), Flat), nothing, Center()),
                                                     ifelse(isa(TZ(), Flat), nothing, Center()))

    ix = interpolator(ii)
    iy = interpolator(jj)
    iz = interpolator(kk)

    i, j, k = (get_node(TX(), Int(ifelse(ix[3] < 0.5, ix[1], ix[2])), grid.Nx),
               get_node(TY(), Int(ifelse(iy[3] < 0.5, iy[1], iy[2])), grid.Ny),
               get_node(TZ(), Int(ifelse(iz[3] < 0.5, iz[1], iz[2])), grid.Nz))

    return i, j, k
end

@inline function apply_tracer_tendency!(::NearestPoint, particles, grid, particle_tendency, tendency, n)
    x = @inbounds particles.x[n]
    y = @inbounds particles.y[n]
    z = @inbounds particles.z[n]

    i, j, k = nearest_node(x, y, z, grid)

    node_volume = volume(i, j, k, grid, Center(), Center(), Center())

    atomic_add!(tendency, i, j, k, particle_tendency / node_volume)

    return nothing
end
using Oceananigans.Fields: AbstractField

# fallback
@inline surface_value(f, args...) = f

struct ContinuousSurfaceFunction{F}
    func :: F
end

@inline function surface_value(f::ContinuousSurfaceFunction, i, j, grid, clock, args...)
    t = clock.time

    x = xnode(i, j, grid.Nz, grid, Center(), Center(), Center())
    y = ynode(i, j, grid.Nz, grid, Center(), Center(), Center())

    return f.func(x, y, t, args...)
end

struct DiscreteSurfaceFuncton{F}
    func :: F
end

@inline surface_value(f::DiscreteSurfaceFuncton, i, j, grid, clock, args...) = f.func(i, j, grid, clock, args...)

# GPU compatible, mainly for fields
@inline surface_value(f::AbstractArray{<:Any, 2}, i, j, grid, clock, args...) = @inbounds f[i, j]
@inline surface_value(f::AbstractArray{<:Any, 3}, i, j, grid, clock, args...) = @inbounds f[i, j, grid.Nz]

# fallback
normalise_surface_function(f; kwargs...) = f

normalise_surface_function(f::Function; discrete_form = false) = discrete_form ? DiscreteSurfaceFuncton(f) : ContinuousSurfaceFunction(f)
using Oceananigans.Fields: AbstractField

# fallback
@inline surface_value(f, args...) = f

struct ContinuousSurfaceFunction{F}
    func :: F
end

@inline Adapt.adapt_structure(to, csf::ContinuousSurfaceFunction) = 
    ContinuousSurfaceFunction(adapt(to, csf.func))

@inline function surface_value(f::ContinuousSurfaceFunction, i, j, grid, clock, args...)
    t = clock.time

    x = xnode(i, j, grid.Nz, grid, Center(), Center(), Center())
    y = ynode(i, j, grid.Nz, grid, Center(), Center(), Center())

    return f.func(x, y, t)
end

struct DiscreteSurfaceFuncton{F}
    func :: F
end

@inline Adapt.adapt_structure(to, csf::DiscreteSurfaceFuncton) = 
    DiscreteSurfaceFuncton(adapt(to, csf.func))

@inline surface_value(f::DiscreteSurfaceFuncton, i, j, grid, clock, args...) = 
    f.func(i, j, grid, clock, args...)

# GPU compatible, mainly for fields
@inline surface_value(f::AbstractArray{<:Any, 2}, i, j, grid, clock, args...) = @inbounds f[i, j]
@inline surface_value(f::AbstractArray{<:Any, 3}, i, j, grid, clock, args...) = @inbounds f[i, j, grid.Nz]

# fallback
normalise_surface_function(f; kwargs...) = f
normalise_surface_function(f::Number; FT, kwargs...) = convert(FT, f)
normalise_surface_function(f::Function; discrete_form = false, kwargs...) = discrete_form ? DiscreteSurfaceFuncton(f) : ContinuousSurfaceFunction(f)

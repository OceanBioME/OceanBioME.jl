using Oceananigans.Fields: fractional_indices, _interpolate, location
using Oceananigans: Center
import Oceananigans.Fields: interpolate

"""
    interpolate(field, x, y, z)
Interpolate `field` to the physical point `(x, y, z)` using trilinear interpolation.
"""
@inline function interpolate(field, x, y, z)
    LX, LY, LZ = location(field)

    if isacolumn(field)
        LZ = Center
        z = field.grid.zᵃᵃᶜ[1]
    end

    i, j, k = fractional_indices(x, y, z, (LX(), LY(), LZ()), field.grid)

    # Convert fractional indices to unit cell coordinates 0 <= (ξ, η, ζ) <=1
    # and integer indices (with 0-based indexing).
    ξ, i = modf(i)
    η, j = modf(j)

    # Convert indices to proper integers and shift to 1-based indexing.
    return _interpolate(field, ξ, η, 0, Int(i+1), Int(j+1), 1)
end

"""
    interpolate(field, LX, LY, LZ, grid, x, y, z)
Interpolate `field` to the physical point `(x, y, z)` using trilinear interpolation. The location of
the field is specified with `(LX, LY, LZ)` and the field is defined on `grid`.
Note that this is a lower-level `interpolate` method defined for use in CPU/GPU kernels.
"""

@inline function interpolate(field, LX, LY, LZ, grid, x, y, z)
    LZ = isacolumn(field) ? Center() : LZ #incase it gets passed LZ=Nothing, e.g. from mean
    i, j, k = fractional_indices(x, y, z, (LX, LY, LZ), grid)

    if isacolumn(field)
        k = 0.0
        ζ = 0.0
    end

    # We use mod and trunc as CUDA.modf is not defined.
    # For why we use Base.unsafe_trunc instead of trunc see:
    # https://github.com/CliMA/Oceananigans.jl/issues/828
    # https://github.com/CliMA/Oceananigans.jl/pull/997
    ξ, i = mod(i, 1), Base.unsafe_trunc(Int, i)
    η, j = mod(j, 1), Base.unsafe_trunc(Int, j)
    ζ, k = mod(k, 1), Base.unsafe_trunc(Int, k)

    return _interpolate(field, ξ, η, ζ, i+1, j+1, k+1)
end
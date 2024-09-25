using KernelAbstractions: @kernel, @index

using Oceananigans.Fields: flatten_node
using Oceananigans.Grids: znode, zspacing

import Oceananigans.Fields: flatten_node

# TODO: move this to Oceananigans
@inline flatten_node(::Nothing, ::Nothing, z) = tuple(z)

@inline shear(z, zₘₓₗ, background_shear, mixed_layer_shear) = ifelse(z <= zₘₓₗ, background_shear, mixed_layer_shear) # Given as 1 in Aumont paper

"""
    ModelLatitude

Returns the latitude specified by the model grid (`y`).
"""
struct ModelLatitude end

"""
    PrescribedLatitude

Returns the prescribed latitude rather than the model grid `y` position.
"""
struct PrescribedLatitude{FT}
    latitude :: FT # °
end

@inline (pl::PrescribedLatitude)(y) = pl.latitude
@inline (::ModelLatitude)(y) = y

"""
    day_length_function(φ, t)

Returns the length of day in seconds at the latitude `φ`, `t`seconds after the start of the year.
"""
@inline function day_length_function(φ, t)
    # as per Forsythe et al., 1995 (https://doi.org/10.1016/0304-3800(94)00034-F)
    p = asind(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (floor(Int, t / day) - 186)))))

    return (24 - 24 / 180 * acosd((sind(0.8333) + sind(φ) * sind(p)) / (cosd(φ) * cosd(p)))) * day / hour
end

"""
    DepthDependantSinkingSpeed(; minimum_speed = 30/day,
                                 maximum_speed = 200/day,
                                 maximum_depth = 500)

Returns sinking speed for particles which sink at `minimum_speed` in the 
surface ocean (the deepest of the mixed and euphotic layers), and accelerate
to `maximum_speed` below that depth and `maximum_depth`.
"""
@kwdef struct DepthDependantSinkingSpeed{FT}
    minimum_speed :: FT = 30/day  # m/s - in NEMO the min and max speeds are both 50m/day
    maximum_speed :: FT = 200/day # m/s
    maximum_depth :: FT = 5000.0  # m
end

# I can't find any explanation as to why this might depend on the euphotic depth
@inline function (p::DepthDependantSinkingSpeed)(i, j, k, grid, mixed_layer_depth, euphotic_depth)
    zₘₓₗ = @inbounds mixed_layer_depth[i, j, k]
    zₑᵤ  = @inbounds euphotic_depth[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())

    w = - p.minimum_speed + (p.maximum_speed - p.minimum_speed) * min(0, z - min(zₘₓₗ, zₑᵤ)) / 5000

    return ifelse(k == grid.Nz + 1, 0, w)
end

# don't actually use this version but might be useful since we can do it
@inline (p::DepthDependantSinkingSpeed)(z, zₘₓₗ, zₑᵤ) = 
    ifelse(z < 0, - p.minimum_speed + (p.maximum_speed - p.minimum_speed) * min(0, z - min(zₘₓₗ, zₑᵤ)) / 5000, 0)

@inline function anoxia_factor(bgc, O₂)
    min_1 = bgc.first_anoxia_threshold
    min_2 = bgc.second_anoxia_threshold

    return min(1, max(0, 0.4 * (min_1 - O₂) / (min_2 + O₂)))
end

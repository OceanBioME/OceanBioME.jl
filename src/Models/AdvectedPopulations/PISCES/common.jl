using Oceananigans.Grids: znode

@inline shear(z, zₘₓₗ, background_shear, mixed_layer_shear) = ifelse(z <= zₘₓₗ, background_shear, mixed_layer_shear) # Given as 1 in Aumont paper

@inline latitude(φ, y) = φ
@inline latitude(::Nothing, y) = y

# we should probably extend this to use DateTime dates at some point
@inline function day_length(φ, t)
    # as per Forsythe et al., 1995 (https://doi.org/10.1016/0304-3800(94)00034-F)
    p = asind(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (floor(Int, t / day) - 186)))))

    return 24 - 24 / 180 * acosd((sind(0.8333) + sind(φ) * sind(p)) / (cosd(φ) * cosd(p)))
end

@kwdef struct DepthDependantSinkingSpeed{FT}
    minimum_speed :: FT = 2/day
    maximum_speed :: FT = 200/day
    maximum_depth :: FT = 5000.0
end

# I can't find any explanation as to why this might depend on the euphotic depth
function (p::DepthDependantSinkingSpeed)(i, j, k, grid, mixed_layer_depth, euphotic_depth)
    zₘₓₗ = @inbounds mixed_layer_depth[i, j, k]
    zₑᵤ  = @inbounds euphotic_depth[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())

    return - p.minimum_speed + (p.maximum_speed - p.minimum_speed) * min(0, z - min(zₘₓₗ, zₑᵤ)) / 5000
end
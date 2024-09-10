@inline shear(z, zₘₓₗ, background_shear, mixed_layer_shear) = ifelse(z <= zₘₓₗ, background_shear, mixed_layer_shear) # Given as 1 in Aumont paper

@inline latitude(φ, y) = φ
@inline latitude(::Nothing, y) = y

# we should probably extend this to use DateTime dates at some point
@inline function day_length(φ, t)
    # as per Forsythe et al., 1995 (https://doi.org/10.1016/0304-3800(94)00034-F)
    p = asind(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (floor(Int, t / day) - 186)))))

    return 24 - 24 / 180 * acosd((sind(0.8333) + sind(φ) * sind(p)) / (cosd(φ) * cosd(p)))
end

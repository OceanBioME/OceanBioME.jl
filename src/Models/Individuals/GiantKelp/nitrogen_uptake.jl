@kwdef struct DroopUptake{CF, FT}
             current_factor :: CF = LowSpeedLimited()

     maximum_ammonia_uptake :: FT = 12 / 0.5 * 24 * 14 / (10^6)# 0.5 from structural_dry_weight_per_area
    ammonia_half_saturation :: FT = 1.3

     maximum_nitrate_uptake :: FT = 10 / 0.5 * 24 * 14 / (10^6)# 0.5 from structural_dry_weight_per_area
    nitrate_half_saturation :: FT = 4.0
end

function potential_ammonia_uptake(uptake::DroopUptake, NH₄, u, v, w)
    jₘ = uptake.maximum_ammonia_uptake
    k = uptake.ammonia_half_saturation

    fₛ = uptake.current_factor(u, v, w)

    return jₘ * fₛ * NH₄ / (NH₄ + k)
end

@inline function (uptake::DroopUptake)(kelp, N, NO₃, NH₄, u, v, w)
    jₘ = uptake.maximum_nitrate_uptake
    k = uptake.nitrate_half_saturation

    Nₘₐ = kelp.maximum_nitrogen_reserve
    Nₘᵢ = kelp.minimum_nitrogen_reserve

    fₛ = uptake.current_factor(u, v, w)

    jNO₃ = max(0, jₘ * fₛ * (Nₘₐ - N) / (Nₘₐ - Nₘᵢ) * NO₃ / (NO₃ + k))

    return jNO₃ # the ammonia "uptake" is actually a substitution of growth consumption
end

@kwdef struct LowSpeedLimited{FT}
    background_factor :: FT = 0.28
      variable_factor :: FT = 0.72
       limiting_speed :: FT = 0.045
end

function (f::LowSpeedLimited)(u, v, w)
    f₁ = f.variable_factor
    f₂ = f.background_factor
    Ū = f.limiting_speed

    U = √(u^2 + v^2 + w^2)
    
    return f₁ * (1 - exp(-U / Ū)) + f₂
end
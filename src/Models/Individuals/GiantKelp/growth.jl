@kwdef struct ReserveAmmoniaGrowth{FT, TL, SL}
            base_growth_rate :: FT = 0.039
maximum_specific_growth_rate :: FT = 0.18
        enhanced_growth_rate :: FT = maximum_specific_growth_rate / (2 * (1 - 0.0126 / 0.0216)) - base_growth_rate
  low_area_enhancement_limit :: FT = 4.5

        temperature_response :: TL = LinearOptimalTemperatureRange()
           seasonal_response :: SL = RateOfDayLengthChange(; latitude = 57.5)
end

function (growth::ReserveAmmoniaGrowth)(kelp, t, A, N, C, T, NH₄, u, v, w)
    μ₀ = growth.base_growth_rate
    μₑ = growth.enhanced_growth_rate
    A₀ = growth.low_area_enhancement_limit

    kₐ = kelp.structural_dry_weight_per_area
    Nₛ = kelp.structural_nitrogen

    Nₘᵢ = kelp.minimum_nitrogen_reserve

    Cₘᵢ = kelp.minimum_carbon_reserve

    μ = (μ₀ + μₑ * exp(-(A / A₀)^2)) # todo verify that Nₘᵢ / Nₘₐ is the correct way up

    fₜ = growth.temperature_response(T)
    fₛ = growth.seasonal_response(t)

    j̃_NH₄ = potential_ammonia_uptake(kelp.nitrogen_uptake, NH₄, u, v, w)
    
    fNH₄ = j̃_NH₄ / kₐ / (N + Nₛ)

    fN = 1 - Nₘᵢ / N
    fC = 1 - Cₘᵢ / C
    
    return μ * fₜ * fₛ * min(fC, max(fN, fNH₄))
end

function carbon_consumption(growth::ReserveAmmoniaGrowth, kelp, t, A, N, C, T, NH₄, u, v, w)
    Cₛ = kelp.structural_carbon

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    return μ * (C + Cₛ)
end

function nitrogen_consumption(growth::ReserveAmmoniaGrowth, kelp, t, A, N, C, T, NH₄, u, v, w)
    Nₛ = kelp.structural_nitrogen

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    consumption = μ * (N + Nₛ)

    jNH₄ = potential_ammonia_uptake(kelp.nitrogen_uptake, NH₄, u, v, w)

    # reserve depletion is substituted by ammonia uptake when possible
    return consumption - min(consumption, jNH₄)
end

function ammonia_replacement(uptake, kelp, N, NH₄, u, v, w, μ)
    Nₛ = kelp.structural_nitrogen

    jNH₄ = potential_ammonia_uptake(uptake, NH₄, u, v, w)

    consumption = μ * (N + Nₛ)

    return min(consumption, jNH₄)
end
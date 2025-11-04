@inline function (kelp::SugarKelp)(::Val{:A}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    ν = erosion(kelp, t, A, N, C, T)

    return A * (μ - ν) / day
end

@inline function (kelp::SugarKelp)(::Val{:N}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    kₐ = kelp.structural_dry_weight_per_area
    Nₛ = kelp.structural_nitrogen

    J = nitrate_uptake(kelp, N, NO₃, u, v, w) +
        ammonia_uptake(kelp, t, A, N, C, T, NH₄, u, v, w)

    e = nitrogen_exudate(kelp, C, T, PAR)

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)
    
    return ((J - e) / kₐ - μ * (N + Nₛ)) / day
end

@inline function (kelp::SugarKelp)(::Val{:C}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    kₐ = kelp.structural_dry_weight_per_area
    Cₛ = kelp.structural_carbon

    P = photosynthesis(kelp, T, PAR)

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    R = respiration(kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)

    e = specific_carbon_exudate(kelp, C)

    return ((P * (1 - e) - R) / kₐ - μ * (C + Cₛ)) / day
end

@inline function growth(kelp, t, A, N, C, T, NH₄, u, v, w)
    f = base_growth_limitation(kelp, t, A, N, C, T)

    # potential ammonia based_growth
    kₐ = kelp.structural_dry_weight_per_area
    Nₛ = kelp.structural_nitrogen

    Nₘ = kelp.minimum_nitrogen_reserve
    Cₘ = kelp.minimum_carbon_reserve
    
    j̃_NH₄ = potential_ammonia_uptake(kelp, NH₄, u, v, w)
    
    μNH₄ = j̃_NH₄ / kₐ / (N + Nₛ)
    
    # internal reserve based growth
    μN = 1 - Nₘ / N
    μC = 1 - Cₘ / C

    return f * min(μC, max(μN, μNH₄))
end

@inline function nitrate_uptake(kelp, N, NO₃, u, v, w)
    J_max_NO₃ = kelp.maximum_nitrate_uptake

    k_NO₃ = kelp.nitrate_half_saturation

    N_max = kelp.maximum_nitrogen_reserve
    N_min = kelp.minimum_nitrogen_reserve

    fᶜ = current_factor(kelp, u, v, w)

    return max(0, J_max_NO₃ * fᶜ * (N_max - N) / (N_max - N_min) * NO₃ / (k_NO₃ + NO₃))
end

@inline function ammonia_uptake(kelp, t, A, N, C, T, NH₄, u, v, w)
    kₐ = kelp.structural_dry_weight_per_area
    Nₛ = kelp.structural_nitrogen

    j̃_NH₄ = potential_ammonia_uptake(kelp, NH₄, u, v, w)

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    return min(j̃_NH₄, μ * kₐ * (N + Nₛ))
end

@inline function potential_ammonia_uptake(kelp, NH₄, u, v, w)
    J_max_NH₄ = kelp.maximum_ammonia_uptake
    k_NH₄ = kelp.ammonia_half_saturation
    
    fᶜ = current_factor(kelp, u, v, w)

    return J_max_NH₄ * fᶜ * NH₄ / (k_NH₄ + NH₄)
end

@inline function photosynthesis(kelp, T, PAR)
    PAR *= day / (3.99e-10 * 545e12) # W / m² / s to einstein / m² / day

    Tk = T + 273.15

    P₁ = kelp.photosynthesis_at_ref_temp_1

    Tₐ  = kelp.photosynthesis_arrhenius_temp
    Tₐₗ = kelp.photosynthesis_low_arrhenius_temp
    Tₐₕ = kelp.photosynthesis_high_arrhenius_temp

    Tₚ  = kelp.photosynthesis_ref_temp_1
    Tₚₗ = kelp.photosynthesis_ref_temp_1
    Tₚₕ = kelp.photosynthesis_high_temp

    α = kelp.photosynthetic_efficiency
    Iₛ = kelp.saturation_irradiance

    maximum_photosynthesis = P₁ * exp(Tₐ / Tₚ - Tₐ / Tk) / (1 + exp(Tₐₗ / Tk - Tₐₗ / Tₚₗ) + exp(Tₐₕ / Tₚₕ - Tₐₕ / Tk))
 
    β = solve_for_light_inhibition(kelp, maximum_photosynthesis)

    pₛ = α * Iₛ / log(1 + α / β)

    return pₛ * (1 - exp(- α * PAR / pₛ)) * exp(-β * PAR / pₛ) 
end


@inline function solve_for_light_inhibition(kelp, Pₘ::FT) where FT
    α = kelp.photosynthetic_efficiency
    Iₛ = kelp.saturation_irradiance

    β₀ = convert(FT, 1e-9)

    return kelp.solver(β_residual, ∂β_maximum_photosynthesis, β₀, (; α, Iₛ, Pₘ))
end

maximum_photosynthesis(α, β) = α / (log(1 + α/β)) * (α / (α + β)) * (β / (α + β)) ^ (β/α)
β_residual(β, p) = (maximum_photosynthesis(p.α, β) - p.Pₘ/p.Iₛ)
∂β_maximum_photosynthesis(β, p) = (p.α * (β / (β + p.α))^(β / p.α) * ((log(p.α / β + 1) * β^2 + p.α * log(p.α / β + 1) * β) * log(β / (β + p.α)) + p.α^2)) / (log(p.α / β + 1)^2 * β * (β + p.α)^2)

@inline function light_inhibition_residual(β, p)
    pₘ = p.maximum_photosynthesis
    α = p.photosynthetic_efficiency
    Iₛ = p.saturation_irradiance

    return pₘ - α * Iₛ / log(1 + α / β) * (α / (α + β)) * (β / (α + β)) ^ (β / α)
end

@inline function respiration(kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)
    Rᵇ = kelp.base_basal_respiration_rate
    Rᵃ = kelp.base_activity_respiration_rate

    Tₐ = kelp.respiration_arrhenius_temp
    T₁ = kelp.respiration_ref_temp_1

    Tk = T + 273.15

    f = exp(Tₐ / T₁ - Tₐ / Tk)

    μₘ = kelp.maximum_specific_growth_rate
    Jₘ = kelp.maximum_nitrate_uptake + kelp.maximum_ammonia_uptake

    J = nitrate_uptake(kelp, N, NO₃, u, v, w) +
        ammonia_uptake(kelp, t, A, N, C, T, NH₄, u, v, w)
    # (basal + activity associated) * temp factor
    return f * (Rᵇ + Rᵃ * (μ / μₘ + J / Jₘ))
end

@inline function specific_carbon_exudate(kelp, C)
    γ = kelp.exudation
    Cₘ = kelp.minimum_carbon_reserve

    return 1 - exp(γ * (Cₘ - C))
end

@inline function nitrogen_exudate(kelp, C, T, PAR)
    CN = kelp.exudation_redfield_ratio

    P = photosynthesis(kelp, T, PAR)
    e = specific_carbon_exudate(kelp, C)

    return P * e * 14 / 12 / CN
end

@inline function erosion(kelp, t, A, N, C, T)
    ε = kelp.erosion_exponent
    ν₀ = kelp.base_erosion_rate

    return ν₀ * exp(ε * A) / (1 + ν₀ * (exp(ε * A) - 1))
end

@inline function current_factor(kelp, u, v, w)
    f₁ = kelp.current_1
    f₂ = kelp.current_2
    f₃ = kelp.current_3

    U = √(u^2 + v^2 + w^2)
    
    return f₁ * (1 - exp(-U / f₃)) + f₂
end

@inline function base_growth_limitation(kelp, t, A, N, C, T)
    fₜ = kelp.temperature_limit(T)
    fₐ = area_limitation(kelp, A)
    fₚ = seasonal_limitation(kelp, t)

    return fₜ * fₐ * fₚ
end

@inline function area_limitation(kelp, A)
    A₀ = kelp.growth_rate_adjustment
    m₁ = kelp.growth_adjustment_1
    m₂ = kelp.growth_adjustment_2

    return m₁ * exp(-(A / A₀)^2) + m₂
end

@kwdef struct LinearOptimalTemperatureRange{FT}
     lower_optimal :: FT = 10.0
     upper_optimal :: FT = 15.0
    lower_gradient :: FT = 1/(lower_optimal + 1.8) #0.08 - I think its important that the minimum cut off (-1.8) is observed, 
    # because the literature doesn't say anything about it growing okay and then suddenly dying off 
    # at low temperature as the origional form suggests, although for higher temperature...
    upper_gradient :: FT = -0.25
end

function (f::LinearOptimalTemperatureRange)(T)
    Tₗ = f.lower_optimal
    Tᵤ = f.upper_optimal

    αₗ = f.lower_gradient
    αᵤ = f.upper_gradient

    return (max(0, αₗ * (T - Tₗ) + 1) * (T < Tₗ) +
            max(0, αᵤ * (T - Tᵤ) + 1) * (T > Tᵤ) +
            one(T) * (Tₗ <= T <= Tᵤ))
end

# the kelp magically know what day of the year it is!!!
@inline function seasonal_limitation(kelp, t)
    φ = kelp.adapted_latitude
    a₁ = kelp.photoperiod_1
    a₂ = kelp.photoperiod_2

    n = floor(Int, mod(t, 364days)/day)

    λ = normed_day_length_change(φ, n)

    return a₁ * (1 + sign(λ) * abs(λ) ^ .5) + a₂
end

# hmmmm
normed_day_length_change(φ, n) = (day_length(φ, n) - day_length(φ, n-1)) / (day_length(φ, 76) - day_length(φ, 75))

@inline function day_length(φ, n)
    n -= 171
    M = mod((356.5291 + 0.98560028 * n), 360)
    C = 1.9148 * sin(M * π / 180) + 0.02 * sin(2 * M * π / 180) + 0.0003 * sin(3 * M * π / 180)
    λ = mod(M + C + 180 + 102.9372, 360)
    δ = asin(sin(λ * π / 180) * sin(23.44 * π / 180))
    ω = (sin(-0.83 * π / 180) * sin(φ * π / 180) * sin(δ)) / (cos(φ * π / 180) * cos(δ))
    return ω / 180
end

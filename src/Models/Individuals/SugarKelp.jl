module SugarKelpModel

using Roots

import OceanBioME.Particles: required_particle_fields, required_tracers, coupled_tracers

# seriously disgsuting number of parameters
@kwdef struct SugarKelp{FT}
    growth_rate_adjustment :: FT = 4.5
    photosynthetic_efficiency :: FT = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60)
    minimum_carbon_reserve :: FT = 0.01
    structural_carbon :: FT = 0.2
    exudation :: FT = 0.5
    erosion :: FT = 0.22
    saturation_irradiance :: FT = 90 * day/ (10 ^ 6)
    structural_dry_weight_per_area :: FT = 0.5
    structural_dry_to_wet_weight :: FT = 0.0785
    carbon_reserve_per_carbon :: FT = 2.1213
    nitrogen_reserve_per_nitrogen :: FT = 2.72
    minimum_nitrogen_reserve :: FT = 0.0126
    maximum_nitrogen_reserve :: FT = 0.0216
    growth_adjustment_2 :: FT = 0.039 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve))
    growth_adjustment_1 :: FT = 0.18 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve)) - growth_adjustment_2
    maximum_specific_growth_rate :: FT = 0.18
    structural_nitrogen :: FT = 0.0146
    photosynthesis_at_ref_temp_1 :: FT = 1.22e-3 * 24
    photosynthesis_at_ref_temp_2 :: FT = 1.3e-3 * 24
    photosynthesis_ref_temp_1 :: FT = 285.0
    photosynthesis_ref_temp_2 :: FT = 288.0
    photoperiod_1 :: FT = 0.85
    photoperiod_2 :: FT = 0.3
    respiration_at_ref_temp_1 :: FT = 2.785e-4 * 24
    respiration_at_ref_temp_2 :: FT = 5.429e-4 * 24
    respiration_ref_temp_1 :: FT = 285.0
    respiration_ref_temp_2 :: FT = 290.0
    photosynthesis_arrhenius_temp :: FT = (1 / photosynthesis_ref_temp_1 - 1 / photosynthesis_ref_temp_2) ^ -1 * log(photosynthesis_at_ref_temp_2 / photosynthesis_at_ref_temp_1)
    photosynthesis_low_temp :: FT = 271.0
    photosynthesis_high_temp :: FT = 296.0
    photosynthesis_high_arrhenius_temp :: FT = 1414.87
    photosynthesis_low_arrhenius_temp :: FT = 4547.89
    respiration_arrhenius_temp :: FT = (1 / respiration_ref_temp_1 - 1 / respiration_ref_temp_2) ^ -1 * log(respiration_at_ref_temp_2 / respiration_at_ref_temp_1)
    current_speed_for_0p65_uptake :: FT = 0.03
    nitrate_half_saturation :: FT = 4.0
    ammonia_half_saturation :: FT = 1.3
    maximum_nitrate_uptake :: FT = 10 * structural_dry_weight_per_area * 24 * 14 / (10^6)
    maximum_ammonia_uptake :: FT = 12 * structural_dry_weight_per_area * 24 * 14 / (10^6)
    current_1 :: FT = 0.72
    current_2 :: FT = 0.28
    current_3 :: FT = 0.045
    base_activity_respiration_rate :: FT = 1.11e-4 * 24
    base_basal_respiration_rate :: FT = 5.57e-5 * 24
    exudation_redfield_ratio :: FT = Inf
end 

@inline required_particle_fields(::SugarKelp) = (:A, :N, :C)
@inline required_tracers(::SugarKelp) = (:u, :v, :w, :NO₃, :NH₄)
@inline coupled_tracers(::SugarKelp) = (:NO₃, :NH₄)

@inline function (kelp::SugarKelp)(::Val{:A}, t, A, N, C, u, v, w, NO₃, NH₄)
    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    ν = errosion(kelp, t, A, N, C, T)

    return A * (μ - ν)
end

@inline function (kelp::SugarKelp)(::Val{:N}, t, A, N, C, u, v, w, NO₃, NH₄)
    kₐ = kelp.structural_dry_weight_per_area
    Nₛ = kelp.structural_nitrogen

    J = nitrate_uptake(kelp, N, NO₃, u, v, w) +
        ammonia_uptake(kelp, t, N, C, T, NH₄, u, v, w)

    e = nitrogen_exudate(kelp, C, T, PAR)

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)
    
    return (J - e) / kₐ - μ * (N + Nₛ)
end

@inline function (kelp::SugarKelp)(::Val{:C}, t, A, N, C, u, v, w, NO₃, NH₄)
    kₐ = kelp.structural_dry_weight_per_area
    Cₛ = kelp.structural_carbon

    P = photosynthesis(kelp, T, PAR)

    R = respiration(kelp, t, N, C, T, NH₄, u, v, w, μ)

    e = specific_carbon_exudate(kelp, C)

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)
    
    return (P * (1 - e) - R) / kₐ - μ * (C + Cₛ)
end

@inline function growth(kelp, t, A, N, C, T, NH₄, u, v, w)
    f = base_growth_limitation(kelp, t, A, N, C, T)

    # potential ammonia based_growth
    kₐ = kelp.structural_dry_weight_per_area
    Nₛ = kelp.structural_nitrogen
    
    j̃_NH₄ = potential_ammonia_uptake(kelp, NH₄, u, v, w)
    
    μNH₄ = j̃_NH₄ / kₐ / (N + Nₛ)
    
    # internal reserve based growth
    μN   = 1 - Nₘ / N
    μC   = 1 - Cₘ / C

    return f * min(μC, max(μN, μNH₄))
end

@inline function nitrate_uptake(kelp, N, NO₃, u, v, w)
    k_NO₃ = kelp.nitrate_half_saturation

    N_max = kelp.maximum_nitrogen_reserve
    N_min = kelp.minimum_nitrogen_reserve

    fᶜ = current_factor(kelp, u, v, w)

    return max(0, J_max_NO₃ * fᶜ * (N_max - N) / (N_max - N_min) * NO₃ / (k_NO₃ + NO₃))
end

@inline function ammonia_uptake(kelp, t, N, C, T, NH₄, u, v, w)
    kₐ = kelp.structural_dry_weight_per_area
    Nₛ = kelp.structural_nitrogen

    j̃_NH₄ = potential_ammonia_uptake(kelp, NH₄, u, v, w)

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    return min(j̃_NH₄, μ * kₐ * (N + Nₛ))
end

@inline function potential_ammonia_uptake(kelp, NH₄, u, v, w)
    fᶜ = current_factor(kelp, u, v, w)

    return J_max_NH₄ * fᶜ * NH₄ / (k_NH₄ + NH₄)
end

@inline function photosynthesis(kelp, T, PAR)
    Tk = T + 273.15

    P₁ = kelp.photosynthesis_at_ref_temp_1

    Tₐ  = kelp.photosynthesis_arrhenius_temp
    Tₐₗ = kelp.photosynthesis_low_arrhenius_temp
    Tₐₕ = kelp.photosynthesis_high_arrhenius_temp

    Tₚ = kelp.photosynthesis_ref_temp_1
    Tₚₗ = kelp.photosynthesis_ref_temp_1
    Tₚₕ = kelp.photosynthesis_high_temp

    maximum_photosynthesis = P₁ * exp(Tₐ / Tₚ - Tₐ / Tk) / (1 + exp(Tₐₗ / Tk - Tₐₗ / Tₚₗ) + exp(Tₐₕ / Tₚₕ - Tₐₕ / Tk))
 
    β = solve_for_light_inhibition(kelp, maximum_photosynthesis)

    pₛ = p.photosynthetic_efficiency * p.saturation_irradiance / log(1 + p.photosynthetic_efficiency / β)

    return pₛ * (1 - exp(- p.photosynthetic_efficiency * PAR / pₛ)) * exp(-β * PAR / pₛ) 
end

# solves `alkalinity_residual` for pH
@inline solve_for_light_inhibition(kelp, maximum_photosynthesis) =
    find_zero(light_inhibition_residual, (0, 0.1), Bisection(); 
              p = (; maximum_photosynthesis, kelp.photosynthetic_efficiency, kelp.saturation_irradiance))

@inline function light_inhibition_residual(β, p)
    pₘ = p.maximum_photosynthesis
    α = p.photosynthetic_efficiency
    Iₛ = p.saturation_irradiance

    return pₘ - α * Iₛ / log(1 + α / β) * (α / (α + β)) * (β / (α + β)) ^ (β / α)
end

@inline function respiration(kelp, t, N, C, T, NH₄, u, v, w, μ)
    Rᵇ = kelp.base_basal_respiration_rate
    Rᵃ = kelp.base_activity_respiration_rate

    Tₐ = kelp.respiration_arrhenius_temp
    T₁ = kelp.respiration_ref_temp_1

    Tk = T + 273.15

    f = exp(Tₐ / T₁ - Tₐ / Tk)

    μₘ = kelp.maximum_specific_growth_rate
    Jₘ = kelp.maximum_nitrate_uptake + kelp.maximum_ammonia_uptake

    J = nitrate_uptake(kelp, N, NO₃, u, v, w) +
        ammonia_uptake(kelp, t, N, C, T, NH₄, u, v, w)
    # (basal + activity associated) * temp factor
    return f * (Rᵇ + Rᵃ * (μ / μₘ + J / Jₘ))
end

@inline function specific_carbon_exudate(kelp, C)
    γ = exudation
    Cₘ = minimum_carbon_reserve

    return 1 - exp(γ * (Cₘ - C))
end

@inline function nitrogen_exudate(kelp, C, T, PAR)
    CN = kelp.exudation_redfield_ratio

    P = photosynthesis(kelp, T, PAR)
    e = specific_carbon_exudate(kelp, C)

    return P * e * 14 / 12 / CN
end
end # module

struct NutrientLimitedProduction{FT}
         base_growth_rate :: FT = 0.6 / day
  temperature_sensetivity :: FT = 1.066
      darkness_tollerance :: FT
initial_slope_of_PI_curve :: FT = 2.0
end

@inline function (μ::NutrientLimitedProduction)(bgc, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR, L)
    bₜ = μ.temperature_sensetivity
    μ₀ = μ.base_growth_rate
    α  = μ.initial_slope_of_PI_curve

    dark_tollerance = μ.dark_tollerance

    θ = IChl / I

    φ = bgc.latitude(y)
    day_length = bgc.day_length(φ, t)

    dark_residence_time = (max(0, zₑᵤ - zₘₓₗ)) ^ 2 / κ

    fₜ = bₜ ^ T

    μᵢ = μ₀ * fₜ

    f₁ = 1.5 * day_length / (day_length + 0.5)

    f₂ = 1 - dark_residence_time / (dark_residence_time + dark_tollerance)

    return μᵢ * f₁ * f₂ * (1 - exp(-α * θ * PAR / (day_length * μᵢ * L))) * L
end

# "new production"
@kwdef struct GrowthRespirationLimitedProduction{FT}
         base_growth_rate :: FT = 0.6 / day
  temperature_sensetivity :: FT = 1.066
          dark_tollerance :: FT
initial_slope_of_PI_curve :: FT = 2.0
   basal_respiration_rate :: FT
    reference_growth_rate :: FT
end

@inline function (μ::GrowthRespirationLimitedProduction)(bgc, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR, L)
    bₜ = μ.temperature_sensetivity
    μ₀ = μ.base_growth_rate
    α  = μ.initial_slope_of_PI_curve
    bᵣ = μ.basal_respiration_rate
    μᵣ = μ.reference_growth_rate

    dark_tollerance = μ.dark_tollerance

    θ = IChl / I

    φ = bgc.latitude(y)
    day_length = bgc.day_length(φ, t)

    dark_residence_time = (max(0, zₑᵤ - zₘₓₗ)) ^ 2 / κ

    fₜ = bₜ ^ T

    μᵢ = μ₀ * fₜ

    f₁ = 1.5 * day_length / (day_length + 0.5)

    f₂ = 1 - dark_residence_time / (dark_residence_time + dark_tollerance)

    return μᵢ * f₁ * f₂ * (1 - exp(-α * θ * PAR / (day_length * (bᵣ + μᵣ) * L))) * L
end

# new method for this if this isn't how you define μ̌
@inline function production_and_energy_assimilation_absorption_ratio(growth_rate, bgc, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR, L)
    α = growth_rate.initial_slope_of_PI_curve

    φ = bgc.latitude(y)
    day_length = bgc.day_length(φ, t)

    f₁ = 1.5 * day_length / (day_length + 0.5)

    μ = growth_rate(bgc, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PARᵢ, L)
    
    μ̌ = μ / f₁

    return μ, 144 * μ̌ * I / (α * IChl * PAR) * day_length
end
    

@inline function base_production_rate(growth_rate, bgc, T)
    bₜ = growth_rate.temperature_sensetivity
    μ₀ = growth_rate.base_growth_rate

    fₜ = bₜ ^ T

    return μ₀ * fₜ
end
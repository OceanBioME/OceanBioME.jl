abstract type BaseProduction end

@inline function (μ::BaseProduction)(phyto, bgc, y, t, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃, L)
    bₜ = μ.temperature_sensetivity
    μ₀ = μ.base_growth_rate
    α₀ = μ.initial_slope_of_PI_curve
    β  = μ.low_light_adaptation

    dark_tollerance = μ.dark_tollerance

    β₁ = phyto.blue_light_absorption
    β₂ = phyto.green_light_absorption
    β₃ = phyto.red_light_absorption

    PAR = β₁ * PAR₁ + β₂ * PAR₂ + β₃ * PAR₃

    φ = bgc.latitude(y)
    day_length = bgc.day_length(φ, t)

    dark_residence_time = max(0, zₑᵤ - zₘₓₗ) ^ 2 / κ

    fₜ = bₜ ^ T

    μᵢ = μ₀ * fₜ

    f₁ = 1.5 * day_length / (day_length + 0.5day)

    f₂ = 1 - dark_residence_time / (dark_residence_time + dark_tollerance)

    α = α₀ * (1 + β * exp(-PAR))

    return μᵢ * f₁ * f₂ * light_limitation(μ, I, IChl, T, PAR, day_length, L, α) * L
end

@inline function (μ::BaseProduction)(phyto, bgc, y, t, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    L, = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    return μ(phyto, bgc, y, t, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃, L)
end

@kwdef struct NutrientLimitedProduction{FT} <: BaseProduction
         base_growth_rate :: FT = 0.6 / day # 1 / s
  temperature_sensetivity :: FT = 1.066     # 
          dark_tollerance :: FT             # s 
initial_slope_of_PI_curve :: FT = 2.0       # 
     low_light_adaptation :: FT = 0.0       # 
end

@inline function light_limitation(μ::NutrientLimitedProduction, I, IChl, T, PAR, day_length, L, α)
    μᵢ = base_production_rate(μ, T)

    θ = IChl / (12 * I + eps(0.0))

    return 1 - exp(-α * θ * PAR / (day_length * μᵢ * L + eps(0.0)))
end

# "new production"
@kwdef struct GrowthRespirationLimitedProduction{FT} <: BaseProduction
         base_growth_rate :: FT = 0.6 / day # 1 / s
  temperature_sensetivity :: FT = 1.066     #
          dark_tollerance :: FT             # s
initial_slope_of_PI_curve :: FT = 2.0       # 
     low_light_adaptation :: FT = 0.0       #
   basal_respiration_rate :: FT = 0.033/day # 1 / s
    reference_growth_rate :: FT = 1.0/day   # 1 / s
end

@inline function light_limitation(μ::GrowthRespirationLimitedProduction, I, IChl, T, PAR, day_length, L, α)
    bᵣ = μ.basal_respiration_rate
    μᵣ = μ.reference_growth_rate

    θ = IChl / (12 * I + eps(0.0))

    return 1 - exp(-α * θ * PAR / (day_length * (bᵣ + μᵣ) * L))
end

# new method for this if this isn't how you define μ̌
@inline function production_and_energy_assimilation_absorption_ratio(growth_rate, phyto, bgc, y, t, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR, PAR₁, PAR₂, PAR₃, L)
    α = growth_rate.initial_slope_of_PI_curve

    φ = bgc.latitude(y)
    day_length = bgc.day_length(φ, t)

    f₁ = 1.5 * day_length / (day_length + 0.5day)

    μ = growth_rate(phyto, bgc, y, t, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃, L)
    μ̌ = μ / f₁

    return μ, 144 * μ̌ * I / (α * IChl * PAR + eps(0.0)) * day_length
end

@inline function base_production_rate(growth_rate, T)
    bₜ = growth_rate.temperature_sensetivity
    μ₀ = growth_rate.base_growth_rate

    fₜ = bₜ ^ T

    return μ₀ * fₜ
end

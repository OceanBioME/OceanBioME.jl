abstract type BaseProduction end

@inline function (μ::BaseProduction)(val_name, i, j, k, grid, bgc, clock, fields, L)
    bₜ = μ.temperature_sensetivity
    μ₀ = μ.base_growth_rate

    dark_tollerance = μ.dark_tollerance

    β₁ = phyto.blue_light_absorption
    β₂ = phyto.green_light_absorption
    β₃ = phyto.red_light_absorption
    
    PAR₁ = @inbounds fields.PAR₁[i, j, k]
    PAR₂ = @inbounds fields.PAR₂[i, j, k]
    PAR₃ = @inbounds fields.PAR₃[i, j, k]

    zₑᵤ = @inbounds fields.zₑᵤ[i, j, k]
    zₘₓₗ = @inbounds fields.zₘₓₗ[i, j, k]

    I, IChl, IFe = phytoplankton_concentration(val_name, i, j, k, fields)

    T = @inbounds fields.T[i, j, k]

    PAR = β₁ * PAR₁ + β₂ * PAR₂ + β₃ * PAR₃

    φ = bgc.latitude(i, j, k, grid)
    day_length = bgc.day_length(φ, t)

    dark_residence_time = max(0, zₑᵤ - zₘₓₗ) ^ 2 / κ

    fₜ = bₜ ^ T

    μᵢ = μ₀ * fₜ

    f₁ = 1.5 * day_length / (day_length + 0.5day)

    f₂ = 1 - dark_residence_time / (dark_residence_time + dark_tollerance)

    α = α₀ * (1 + β * exp(-PAR))

    fₗ = light_limitation(μ, I, IChl, T, PAR, day_length, L, α)

    return μᵢ * f₁ * f₂ * fₗ * L
end

"""
    NutrientLimitedProduction

`BaseProduction` with light limitation moderated by nutrient availability. This is
the "origional" PISCES phytoplankton growth rate model. Growth rate is of the form:

```math
μ = μ₁f₁(τᵈ)f₂(zₘₓₗ)(1-exp(-α θᶜʰˡ PAR / τ μ₀ L)) L.
```

Keyword Arguments
=================
- `base_growth_rate`: the base growth rate, μ₀, in (1/s)
- `temperatrue_sensetivity`: temperature sensetivity parameter, b, giving μ₁ = μ₀ bᵀ where T is temperature
- `dark_tollerance`: the time that the phytoplankton survives in darkness below the euphotic layer, τᵈ (s)
- `initial_slope_of_PI_curve`: the relationship between photosynthesis and irradiance, α₀ (1/W/m²)
- `low_light_adaptation`: factor increasing the sensetivity of photosynthesis to irradiance, β,
   giving α = α₀(1 + exp(-PAR)), typically set to zero

"""
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

"""
    NutrientLimitedProduction

`BaseProduction` with light limitation moderated by nutrient availability. This is
the "new production" PISCES phytoplankton growth rate model. Growth rate is of the form:

```math
μ = μ₁f₁(τ)f₂(zₘₓₗ)(1-exp(-α θᶜʰˡ PAR / τ (bᵣ + μᵣ))) L.
```

Keyword Arguments
=================
- `base_growth_rate`: the base growth rate, μ₀, in (1/s)
- `temperatrue_sensetivity`: temperature sensetivity parameter, b, giving μ₁ = μ₀ bᵀ where T is temperature
- `dark_tollerance`: the time that the phytoplankton survives in darkness below the euphotic layer, τᵈ (s)
- `initial_slope_of_PI_curve`: the relationship between photosynthesis and irradiance, α₀ (1/W/m²)
- `low_light_adaptation`: factor increasing the sensetivity of photosynthesis to irradiance, β,
   giving α = α₀(1 + exp(-PAR)), typically set to zero
- `basal_respiration_rate`: reference respiration rate, bᵣ (1/s)
- `reference_growth_rate`: reference growth rate, μᵣ (1/s)
   
"""
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

    return 1 - exp(-α * θ * PAR / (day_length * (bᵣ + μᵣ)))
end

@inline function production_and_energy_assimilation_absorption_ratio(phyto, val_name, bgc, i, j, k, grid, bgc, clock, fields)
    α₀ = growth_rate.initial_slope_of_PI_curve
    β  = growth_rate.low_light_adaptation
    β₁ = phyto.blue_light_absorption
    β₂ = phyto.green_light_absorption
    β₃ = phyto.red_light_absorption
    
    I, IChl, IFe = phytoplankton_concentration(val_name, i, j, k, fields)
    
    PAR₁ = @inbounds fields.PAR₁[i, j, k]
    PAR₂ = @inbounds fields.PAR₂[i, j, k]
    PAR₃ = @inbounds fields.PAR₃[i, j, k]

    PAR = β₁ * PAR₁ + β₂ * PAR₂ + β₃ * PAR₃

    φ = bgc.latitude(i, j, k, grid)
    day_length = bgc.day_length(φ, t)

    f₁ = 1.5 * day_length / (day_length + 0.5day)

    μ = growth_rate(val_name, i, j, k, grid, bgc, clock, fields, L)

    μ̌ = μ / f₁ * day_length 

    α = α₀ * (1 + β * exp(-PAR))

    return μ, 12 * μ̌ * I / (α * IChl * PAR + eps(0.0)) * L # (1 / s, unitless)
end

@inline function base_production_rate(growth_rate, T)
    bₜ = growth_rate.temperature_sensetivity
    μ₀ = growth_rate.base_growth_rate

    fₜ = bₜ ^ T

    return μ₀ * fₜ
end

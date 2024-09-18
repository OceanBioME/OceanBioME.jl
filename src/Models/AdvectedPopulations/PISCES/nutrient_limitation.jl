@kwdef struct NitrogenIronPhosphateSilicateLimitation{FT, BT}
    minimum_ammonium_half_saturation :: FT
     minimum_nitrate_half_saturation :: FT
   minimum_phosphate_half_saturation :: FT
       threshold_for_size_dependency :: FT = 1.0
                          size_ratio :: FT = 3.0
                  optimal_iron_quota :: FT = 7.0e-3
                    silicate_limited :: BT
    minimum_silicate_half_saturation :: FT = 1.0
  silicate_half_saturation_parameter :: FT = 16.6
     half_saturation_for_iron_uptake :: FT
end

@inline function size_factor(L, I)
    Iₘ  = L.threshold_for_size_dependency
    S   = L.size_ratio

    I₁ = min(I, Iₘ)
    I₂ = max(0, I - Iₘ)

    return (I₁ + S * I₂) / (I₁ + I₂ + eps(0.0))
end

@inline function (L::NitrogenIronPhosphateSilicateLimitation)(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)
    kₙₒ = L.minimum_nitrate_half_saturation
    kₙₕ = L.minimum_ammonium_half_saturation
    kₚ  = L.minimum_phosphate_half_saturation
    kₛᵢ = L.minimum_silicate_half_saturation
    pk  = L.silicate_half_saturation_parameter

    θₒ  = L.optimal_iron_quota

    # quotas
    θFe  = ifelse(I == 0, 0, IFe / (I + eps(0.0)))
    θChl = ifelse(I == 0, 0, IChl / (12 * I + eps(0.0)))

    K̄ = size_factor(L, I)

    Kₙₒ = kₙₒ * K̄
    Kₙₕ = kₙₕ * K̄
    Kₚ  = kₚ  * K̄

    # nitrogen limitation
    LNO₃ = nitrogen_limitation(NO₃, NH₄, Kₙₒ, Kₙₕ)
    LNH₄ = nitrogen_limitation(NH₄, NO₃, Kₙₕ, Kₙₒ)

    LN = LNO₃ + LNH₄

    # phosphate limitation
    LPO₄ = PO₄ / (PO₄ + Kₚ + eps(0.0))

    # iron limitation
    # Flynn and Hipkin (1999) - photosphotosyntheis, respiration (?), nitrate reduction 
    θₘ = 0.0016 / 55.85 * θChl + 1.5 * 1.21e-5 * 14 / (55.85 * 7.625) * LN + 1.15e-4 * 14 / (55.85 * 7.625) * LNO₃

    LFe = min(1, max(0, (θFe - θₘ) / θₒ))

    # silicate limitation
    KSi = kₛᵢ + 7 * Si′^2 / (pk^2 + Si′^2)
    LSi = Si / (Si + KSi)
    LSi = ifelse(L.silicate_limited, LSi, Inf)

    # don't always need the other arguments but they can be got like L, = ... or _, LFe = ..
    return min(LN, LPO₄, LFe, LSi), LFe, LPO₄, LN, LNO₃, LNH₄
end

@inline nitrogen_limitation(N₁, N₂, K₁, K₂) = (K₂ * N₁) / (K₁ * K₂ + K₁ * N₂ + K₂ * N₁ + eps(0.0))

@inline function iron_uptake_limitation(L, I, Fe)
    k = L.half_saturation_for_iron_uptake

    K = k * size_factor(L, I)

    return Fe / (Fe + K + eps(0.0))
end

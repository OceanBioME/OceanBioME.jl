"""
    NitrogenIronPhosphateSilicateLimitation

Holds the parameters for growth limitation by nitrogen (NO₃ and NH₄),
iron (Fe), phosphate PO₄, and (optionally) silicate (Si) availability.

Silicate limitation may be turned off (e.g. for nanophytoplankton) by
setting `silicate_limited=false`.
"""
@kwdef struct NitrogenIronPhosphateSilicateLimitation{FT}
    minimum_ammonium_half_saturation :: FT          # mmol N / m³
     minimum_nitrate_half_saturation :: FT          # mmol N / m³
   minimum_phosphate_half_saturation :: FT          # mmol P / m³
                  optimal_iron_quota :: FT = 0.007  # μmol Fe / mmol C
                    silicate_limited :: Bool        # Bool
    minimum_silicate_half_saturation :: FT = 1.0    # mmol Si / m³
  silicate_half_saturation_parameter :: FT = 16.6   # mmol Si / m³
end

@inline function (L::NitrogenIronPhosphateSilicateLimitation)(phyto, I, IChl, IFe, NO₃, NH₄, PO₄, Si, Si′)
    kₙₒ = L.minimum_nitrate_half_saturation
    kₙₕ = L.minimum_ammonium_half_saturation
    kₚ  = L.minimum_phosphate_half_saturation
    kₛᵢ = L.minimum_silicate_half_saturation
    pk  = L.silicate_half_saturation_parameter

    θₒ  = L.optimal_iron_quota

    θFe  = ifelse(I == 0, 0, IFe / (I + eps(0.0)))
    θChl = ifelse(I == 0, 0, IChl / (12 * I + eps(0.0)))

    K̄ = size_factor(phyto, I)

    Kₙₒ = kₙₒ * K̄
    Kₙₕ = kₙₕ * K̄
    Kₚ  = kₚ  * K̄
    Kₛᵢ = kₛᵢ * K̄

    LNO₃ = nitrogen_limitation(NO₃, NH₄, Kₙₒ, Kₙₕ)
    LNH₄ = nitrogen_limitation(NH₄, NO₃, Kₙₕ, Kₙₒ)

    LN = LNO₃ + LNH₄

    LPO₄ = PO₄ / (PO₄ + Kₚ + eps(0.0))

    θₘ = 10^3 * (0.0016 / 55.85 * 12 * θChl + 1.5 * 1.21e-5 * 14 / (55.85 * 7.625) * LN + 1.15e-4 * 14 / (55.85 * 7.625) * LNO₃)

    LFe = min(1, max(0, (θFe - θₘ) / θₒ))

    KSi = Kₛᵢ + 7 * Si′^2 / (pk^2 + Si′^2)
    LSi = Si / (Si + KSi)
    LSi = ifelse(L.silicate_limited, LSi, Inf)

    return min(LN, LPO₄, LFe, LSi), LFe, LPO₄, LN, LNO₃, LNH₄
end

@inline function (L::NitrogenIronPhosphateSilicateLimitation)(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields)
    kₙₒ = L.minimum_nitrate_half_saturation
    kₙₕ = L.minimum_ammonium_half_saturation
    kₚ  = L.minimum_phosphate_half_saturation
    kₛᵢ = L.minimum_silicate_half_saturation
    pk  = L.silicate_half_saturation_parameter

    θₒ  = L.optimal_iron_quota

    I, IChl, IFe = phytoplankton_concentrations(val_name, i, j, k, fields)

    NO₃ = @inbounds fields.NO₃[i, j, k]
    NH₄ = @inbounds fields.NH₄[i, j, k]
    PO₄ = @inbounds fields.PO₄[i, j, k]
    Si  = @inbounds  fields.Si[i, j, k]

    Si′ = @inbounds bgc.silicate_climatology[i, j, k]

    # quotas
    θFe  = ifelse(I == 0, 0, IFe / (I + eps(0.0)))
    θChl = ifelse(I == 0, 0, IChl / (12 * I + eps(0.0)))

    K̄ = size_factor(phyto, I)

    Kₙₒ = kₙₒ * K̄
    Kₙₕ = kₙₕ * K̄
    Kₚ  = kₚ  * K̄
    Kₛᵢ = kₛᵢ * K̄

    # nitrogen limitation
    LNO₃ = nitrogen_limitation(NO₃, NH₄, Kₙₒ, Kₙₕ)
    LNH₄ = nitrogen_limitation(NH₄, NO₃, Kₙₕ, Kₙₒ)

    LN = LNO₃ + LNH₄

    # phosphate limitation
    LPO₄ = PO₄ / (PO₄ + Kₚ + eps(0.0))

    # iron limitation
    # Flynn and Hipkin (1999) - photosynthesis, respiration (?), nitrate reduction 
    θₘ = 10^3 * (0.0016 / 55.85 * 12 * θChl + 1.5 * 1.21e-5 * 14 / (55.85 * 7.625) * LN + 1.15e-4 * 14 / (55.85 * 7.625) * LNO₃) # 1 / 1 to 1/10^3 / 1
    
    LFe = min(1, max(0, (θFe - θₘ) / θₒ))

    # silicate limitation
    KSi = Kₛᵢ + 7 * Si′^2 / (pk^2 + Si′^2)
    LSi = Si / (Si + KSi)
    LSi = ifelse(L.silicate_limited, LSi, Inf)

    # don't always need the other arguments but they can be got like L, = ... or _, LFe = ..
    return min(LN, LPO₄, LFe, LSi), LFe, LPO₄, LN, LNO₃, LNH₄
end

@inline nitrogen_limitation(N₁, N₂, K₁, K₂) = (K₂ * N₁) / (K₁ * K₂ + K₁ * N₂ + K₂ * N₁ + eps(0.0))

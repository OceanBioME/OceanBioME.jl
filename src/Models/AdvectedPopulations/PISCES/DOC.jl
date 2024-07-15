@inline function Rᵤₚ(M, T)
    σᴹ = bgc.non_assimilated_fraction[2]
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton[2]
    mᴹ = bgc.phytoplankton_mortality_rate[2]
    bₘ = bgc.temperature_sensitivity_term[2]
    return (1 - σᴹ - eₘₐₓᴹ)*(1)/(1-eₘₐₓᴹ)*mᴹ*bₘ^T*M^2  #30b
end

@inline function Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, Bactᵣₑ)
    O₂ᵘᵗ = bgc.OC_for_ammonium_based_processes
    λ_DOC = bgc.remineralisation_rate_of_DOC
    bₚ = bgc.temperature_sensitivity_of_growth

    Lₗᵢₘᵇᵃᶜᵗ = Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe)[2]

    return min(O₂/O₂ᵘᵗ, λ_DOC*bₚ^T(1 - ΔO₂(O₂)) * Lₗᵢₘᵇᵃᶜᵗ * (Bact)/(Bactᵣₑ) * DOC) #33a
end

@inline function Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, Bactᵣₑ)
    λ_DOC = bgc.remineralisation_rate_of_DOC
    rₙₒ₃¹ = bgc.CN_ratio_of_denitrification
    bₚ = bgc.temperature_sensitivity_of_growth

    Lₗᵢₘᵇᵃᶜᵗ = Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe)[2]

    return min(NO₃/rₙₒ₃¹, λ_DOC*bₚ^T* ΔO₂(O₂)* Lₗᵢₘᵇᵃᶜᵗ*(Bact)/(Bactᵣₑ) * DOC) #33b
end

@inline function Bact()

        # not sure how to do this
    return
end

@inline function Φ(DOC, POC, GOC, sh)
    a₁ = bgc.aggregation_rate_of_DOC_to_POC_1
    a₂ = bgc.aggregation_rate_of_DOC_to_POC_2
    a₃ = bgc.aggregation_rate_of_DOC_to_GOC_3
    a₄ = bgc.aggregation_rate_of_DOC_to_POC_4
    a₅ = bgc.aggregation_rate_of_DOC_to_POC_5

    Φ₁ᴰᴼᶜ = sh * (a₁*DOC + a₂*POC)*DOC  #36a
    Φ₂ᴰᴼᶜ = sh * (a₃*GOC) * DOC         #36b
    Φ₃ᴰᴼᶜ = (a₄*POC + a₅*DOC)*DOC       #36c

    return Φ₁ᴰᴼᶜ, Φ₂ᴰᴼᶜ, Φ₃ᴰᴼᶜ
end

@inline function Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe)
   
    Kₚₒ₄ᵇᵃᶜᵗ = bgc.PO4_half_saturation_const_for_DOC_remin
    Kₙₒ₃ᵇᵃᶜᵗ = bgc.NO3_half_saturation_const_for_DOC_remin
    Kₙₕ₄ᵇᵃᶜᵗ = bgc.NH4_half_saturation_const_for_DOC_remin
    K_Feᵇᵃᶜᵗ = bgc.Fe_half_saturation_const_for_DOC_remin
    K_DOC = bgc.half_saturation_const_for_DOC_remin

    L_DOCᵇᵃᶜᵗ = K_mondo(DOC, K_DOC) #34b
    L_Feᵇᵃᶜᵗ = K_mondo(bFe, K_Feᵇᵃᶜᵗ) #34d
    Lₚₒ₄ᵇᵃᶜᵗ = K_mondo(PO₄, Kₚₒ₄ᵇᵃᶜᵗ) #34e

    Lₙₕ₄ᵇᵃᶜᵗ = L_NH₄(NO₃, NH₄, Kₙₒ₃ᵇᵃᶜᵗ, Kₙₕ₄ᵇᵃᶜᵗ) #34g
    Lₙₒ₃ᵇᵃᶜᵗ = L_NO₃(NO₃, NH₄, Kₙₒ₃ᵇᵃᶜᵗ, Kₙₕ₄ᵇᵃᶜᵗ) #34h
    Lₙ = Lₙₒ₃ᵇᵃᶜᵗ + Lₙₕ₄ᵇᵃᶜᵗ         #34f

    Lₗᵢₘᵇᵃᶜᵗ = min(Lₙₕ₄ᵇᵃᶜᵗ, Lₚₒ₄ᵇᵃᶜᵗ, L_Feᵇᵃᶜᵗ) #34c
    Lᵇᵃᶜᵗ = Lₗᵢₘᵇᵃᶜᵗ*L_DOCᵇᵃᶜᵗ #34a

    return Lᵇᵃᶜᵗ, Lₗᵢₘᵇᵃᶜᵗ
end


@inline function (pisces::PISCES)(::Val{:DOC}, x, y, z, t, P, D, Pᶜʰˡ, Dᶜʰˡ, N, Fe, O₂, NO₃, PARᴾ, PARᴰ, Z, M, POC, GOC, T, L_day, zₘₓₗ, zₑᵤ)
    [γᶻ, γᴹ] = bgc.excretion_as_DOM
    [σᶻ, σᴹ] = bgc.non_assimilated_fraction
    [δᴾ, δᴰ] = bgc.exudation_of_DOC
    [eₘₐₓᶻ, eₘₐₓᴹ] = bgc.max_growth_efficiency_of_zooplankton
    [αᴾ, αᴰ] = bgc.initial_slope_of_PI_curve

    g_FF = bgc.flux_feeding_rate
    w_GOCᵐⁱⁿ = bgc.min_sinking_speed_of_GOC
    bₘ = bgc.temperature_sensitivity_term[2]

    ∑ᵢgᵢᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ  = grazingᶻ(P, D, POC, T) 
    ∑ᵢgᵢᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ = grazingᴹ(P, D, Z, POC, T)

    zₘₐₓ = max(zₑᵤ, zₘₓₗ)   #41a
    w_GOC = w_GOCᵐⁱⁿ + (200 - w_GOCᵐⁱⁿ)*(max(0, z-zₘₐₓ))/(5000) #41b
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC #29b

    t_darkᴾ =
    t_darkᴰ =

    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[1]
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[1]

    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ)
    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, N, Fe, P, D, POC, 1, Z)
    eᴹ = eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ, N, Fe, P, D, POC, Z, M)

    λₚₒ¹ = 
    Rᵤₚᴹ = Rᵤₚ(M, T)

    zₘₐₓ = max(zₑᵤ, zₘₓₗ) #35a
    Bact = Bact()
    Bactᵣₑ
  
    Remin = Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, Bactᵣₑ)
    Denit = Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, Bactᵣₑ)

    Φ₁ᴰᴼᶜ, Φ₂ᴰᴼᶜ, Φ₃ᴰᴼᶜ = Φ(DOC, POC, GOC, sh)

    return (1 - γᶻ)*(1 - eᶻ - σᶻ)*∑ᵢgᵢᶻ*Z + (1 - γᴹ)*(1 - eᴹ - σᴹ)*(∑ᵢgᵢᴹ + g_GOC_FFᴹ)*M + δᴰ*μᴰ*D + δᴾ*μᴾ*P + λₚₒ¹*POC + (1 - γᴹ)*Rᵤₚᴹ - Remin - Denit - Φ₁ᴰᴼᶜ - Φ₂ᴰᴼᶜ - Φ₃ᴰᴼᶜ #32
end
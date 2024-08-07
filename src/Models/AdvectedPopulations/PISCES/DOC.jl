@inline function Rᵤₚ(M, T, bgc) #third term has small magnitude, as mᴹ per day
    σᴹ = bgc.non_assimilated_fraction.M
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    mᴹ = bgc.zooplankton_quadratic_mortality.M
    bₘ = bgc.temperature_sensitivity_term.M
    return (1 - σᴹ)*mᴹ*(bₘ^T)*M^2  #30b
end

@inline function Pᵤₚ(M, T, bgc)
    σᴹ = bgc.non_assimilated_fraction.M
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    mᴹ = bgc.zooplankton_quadratic_mortality.M
    bₘ = bgc.temperature_sensitivity_term.M
    return σᴹ*mᴹ*(bₘ^T)*M^2      #30a
end


@inline function get_Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc)
    O₂ᵘᵗ = bgc.OC_for_ammonium_based_processes
    λ_DOC = bgc.remineralisation_rate_of_DOC
    bₚ = bgc.temperature_sensitivity_of_growth
    Bactᵣₑ = bgc.bacterial_reference

    Lₗᵢₘᵇᵃᶜᵗ = Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe, bgc)[2]

    #min((O₂)/(O₂ᵘᵗ), λ_DOC*bₚ^T*(1 - ΔO₂(O₂, bgc)) * Lₗᵢₘᵇᵃᶜᵗ * (Bact)/(Bactᵣₑ) * DOC), definition did not make sense with dimensions
    return λ_DOC*bₚ^T*(1 - ΔO₂(O₂, bgc)) * Lₗᵢₘᵇᵃᶜᵗ * (Bact)/(Bactᵣₑ + eps(0.0)) * DOC #33a
end

@inline function get_Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc)
    λ_DOC = bgc.remineralisation_rate_of_DOC
    rₙₒ₃¹ = bgc.CN_ratio_of_denitrification
    bₚ = bgc.temperature_sensitivity_of_growth
    Bactᵣₑ = bgc.bacterial_reference

    Lₗᵢₘᵇᵃᶜᵗ = Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe, bgc)[2]

    #min(NO₃/rₙₒ₃¹, λ_DOC*bₚ^T* ΔO₂(O₂, bgc)* Lₗᵢₘᵇᵃᶜᵗ*(Bact)/(Bactᵣₑ) * DOC), definition did not make sense with dimensions
    return λ_DOC*bₚ^T* ΔO₂(O₂, bgc)* Lₗᵢₘᵇᵃᶜᵗ*(Bact)/(Bactᵣₑ + eps(0.0)) * DOC #33b
end

@inline get_Bact(zₘₐₓ, z, Z, M) = ifelse(z <= zₘₐₓ, min(0.7*(Z + 2*M), 4), min(0.7*(Z + 2*M), 4)*(zₘₐₓ/(z + eps(0.0))^0.683))  #35b

@inline function Φᴰᴼᶜ(DOC, POC, GOC, sh, bgc)
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

@inline function λ¹(T, O₂, bgc) #has small magnitude as λₚₒ per day
    λₚₒ= bgc.degradation_rate_of_POC
    bₚ = bgc.temperature_sensitivity_of_growth
    return λₚₒ*bₚ^T*(1 - 0.45*ΔO₂(O₂, bgc))  #38
end

@inline function Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe, bgc)
   
    Kₚₒ₄ᵇᵃᶜᵗ = bgc.PO4_half_saturation_const_for_DOC_remin
    Kₙₒ₃ᵇᵃᶜᵗ = bgc.NO3_half_saturation_const_for_DOC_remin
    Kₙₕ₄ᵇᵃᶜᵗ = bgc.NH4_half_saturation_const_for_DOC_remin
    K_Feᵇᵃᶜᵗ = bgc.Fe_half_saturation_const_for_DOC_remin
    K_DOC = bgc.half_saturation_const_for_DOC_remin

    L_DOCᵇᵃᶜᵗ = concentration_limitation(DOC, K_DOC) #34b
    L_Feᵇᵃᶜᵗ = concentration_limitation(bFe, K_Feᵇᵃᶜᵗ) #34d
    Lₚₒ₄ᵇᵃᶜᵗ = concentration_limitation(PO₄, Kₚₒ₄ᵇᵃᶜᵗ) #34e

    Lₙₕ₄ᵇᵃᶜᵗ = L_NH₄(NO₃, NH₄, Kₙₒ₃ᵇᵃᶜᵗ, Kₙₕ₄ᵇᵃᶜᵗ) #34g
    Lₙₒ₃ᵇᵃᶜᵗ = L_NO₃(NO₃, NH₄, Kₙₒ₃ᵇᵃᶜᵗ, Kₙₕ₄ᵇᵃᶜᵗ) #34h
    Lₙᵇᵃᶜᵗ = Lₙₒ₃ᵇᵃᶜᵗ + Lₙₕ₄ᵇᵃᶜᵗ         #34f
    #Lₙᵇᵃᶜᵗ is not used...
    Lₗᵢₘᵇᵃᶜᵗ = min(Lₙₕ₄ᵇᵃᶜᵗ, Lₚₒ₄ᵇᵃᶜᵗ, L_Feᵇᵃᶜᵗ) #34c
    Lᵇᵃᶜᵗᵣ = Lₗᵢₘᵇᵃᶜᵗ*L_DOCᵇᵃᶜᵗ #34a        #Lᵇᵃᶜᵗᵣ to avoid method call error

    return Lᵇᵃᶜᵗᵣ, Lₗᵢₘᵇᵃᶜᵗ
end


@inline function (bgc::PISCES)(::Val{:DOC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³)
    γᶻ = bgc.excretion_as_DOM.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᶻ = bgc.non_assimilated_fraction.Z
    σᴹ = bgc.non_assimilated_fraction.M
    δᴾ = bgc.exudation_of_DOC.P
    δᴰ = bgc.exudation_of_DOC.D
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton.Z
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    αᴾ= bgc.initial_slope_of_PI_curve.P
    αᴰ = bgc.initial_slope_of_PI_curve.D
    wₚₒ = bgc.sinking_speed_of_POC

    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)


    g_FF = bgc.flux_feeding_rate
    w_GOCᵐⁱⁿ = bgc.min_sinking_speed_of_GOC
    bₘ = bgc.temperature_sensitivity_term.M

    ∑ᵢgᵢᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ  = get_grazingᶻ(P, D, POC, T, bgc) 
    ∑ᵢgᵢᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)

    zₘₐₓ = max(zₑᵤ, zₘₓₗ)   #41a
    w_GOC = w_GOCᵐⁱⁿ + (200 - w_GOCᵐⁱⁿ)*(max(0, z-zₘₐₓ))/(5000) #41b
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC #29b
    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC

    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = get_PARᴾ(PAR¹, PAR², PAR³, bgc)
    PARᴰ = get_PARᴰ(PAR¹, PAR², PAR³, bgc)

    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]

    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc)
    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eᴹ = eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    λₚₒ¹ = λ¹(T, O₂, bgc)
    Rᵤₚᴹ = Rᵤₚ(M, T, bgc)

    Bact = get_Bact(zₘₐₓ, z, Z, M)

    bFe = Fe #defined in previous PISCES model
    sh = get_sh(z, zₘₓₗ)
  
    Remin = get_Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc)
    Denit = get_Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc)

    Φ₁ᴰᴼᶜ, Φ₂ᴰᴼᶜ, Φ₃ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh, bgc)

    return ((1 - γᶻ)*(1 - eᶻ - σᶻ)*∑ᵢgᵢᶻ*Z + (1 - γᴹ)*(1 - eᴹ - σᴹ)*(∑ᵢgᵢᴹ + g_GOC_FFᴹ + gₚₒ_FFᴹ)*M + 
           δᴰ*μᴰ*D + δᴾ*μᴾ*P + λₚₒ¹*POC + (1 - γᴹ)*Rᵤₚᴹ - Remin - Denit - Φ₁ᴰᴼᶜ - Φ₂ᴰᴼᶜ - Φ₃ᴰᴼᶜ) #32
end #changed this to include gₚₒ_FF
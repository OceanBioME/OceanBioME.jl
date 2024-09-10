# This document contains functions for:
    # grazing
    # gross growth efficiency
    # Z and M forcing

    # Checked all eqns
    # Simplifications possible
    # Could simplify eₙᴶ functions

#Mesozooplankton are grazed by upper trophic levels. Carbon is returned to the system through fecal pellets and respiration.
#Respiration and excretion from upper trophic levels.
@inline function upper_respiration(M, T, bgc) #third term has small magnitude, as mᴹ per day
    σᴹ = bgc.non_assimilated_fraction.M
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    mᴹ = bgc.zooplankton_quadratic_mortality.M
    bₘ = bgc.temperature_sensitivity_term.M
    return (1 - σᴹ - eₘₐₓᴹ)*(1/(1-eₘₐₓᴹ))*mᴹ*(bₘ^T)*M^2  #30b
end

#Fecal pellets from upper trophic levels
@inline function production_of_fecal_pellets(M, T, bgc)
    σᴹ = bgc.non_assimilated_fraction.M
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    mᴹ = bgc.zooplankton_quadratic_mortality.M
    bₘ = bgc.temperature_sensitivity_term.M
    return σᴹ*mᴹ*(1/(1-eₘₐₓᴹ))*(bₘ^T)*M^2      #30a
end

#Zooplankton graze on P, D, POC. We return a list of grazing of Z on each prey and a sum of grazing terms.
@inline function grazing_Z(P, D, POC, T, bgc) 
    pₚᶻ = bgc.preference_for_nanophytoplankton.Z
    p_Dᶻ = bgc.preference_for_diatoms.Z
    pₚₒᶻ = bgc.preference_for_POC.Z
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Fₜₕᵣₑₛₕᶻ = bgc.food_threshold_for_zooplankton.Z
    gₘₐₓᶻ = bgc.max_grazing_rate.Z
    K_Gᶻ = bgc.half_saturation_const_for_grazing.Z
    b_Z = bgc.temperature_sensitivity_term.Z

    F = pₚᶻ*max(0, P - Jₜₕᵣₑₛₕᶻ) + p_Dᶻ*max(0, D - Jₜₕᵣₑₛₕᶻ) + pₚₒᶻ*max(0, POC - Jₜₕᵣₑₛₕᶻ)
    Fₗᵢₘ = max(0, F - min(0.5*F, Fₜₕᵣₑₛₕᶻ))

    grazing_arg = gₘₐₓᶻ*(b_Z^T)*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᶻ + pₚᶻ*P + p_Dᶻ*D + pₚₒᶻ*POC + eps(0.0)))

    gₚᶻ = (pₚᶻ*max(0, P - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    g_Dᶻ = (p_Dᶻ*max(0, D - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    gₚₒᶻ = (pₚₒᶻ*max(0, POC - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    ∑gᶻ= gₚᶻ + g_Dᶻ + gₚₒᶻ  #Sum grazing rates on each prey species for microzooplankton

    return ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ #eq 26a
end

#Mesozooplankton graze on P, D, POC, Z. We return a list of grazing of M on each prey and a sum of grazing terms.
@inline function grazing_M(P, D, Z, POC, T, bgc) #eq 26a
    pₚᴹ = bgc.preference_for_nanophytoplankton.M
    p_Dᴹ = bgc.preference_for_diatoms.M
    pₚₒᴹ = bgc.preference_for_POC.M
    p_Zᴹ = bgc.preference_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    Fₜₕᵣₑₛₕᴹ = bgc.food_threshold_for_zooplankton.M
    gₘₐₓᴹ = bgc.max_grazing_rate.M
    K_Gᴹ = bgc.half_saturation_const_for_grazing.M
    bₘ = bgc.temperature_sensitivity_term.M
    
    F = pₚᴹ*max(0, P - Jₜₕᵣₑₛₕᴹ) + p_Dᴹ*max(0, D - Jₜₕᵣₑₛₕᴹ) + pₚₒᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ) + p_Zᴹ*max(0, Z - Jₜₕᵣₑₛₕᴹ)
    Fₗᵢₘ = max(0, F - min(0.5*F, Fₜₕᵣₑₛₕᴹ))

    grazing_arg =  gₘₐₓᴹ*(bₘ^T)*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᴹ + pₚᴹ*P + p_Dᴹ*D + pₚₒᴹ*POC + p_Zᴹ*Z + eps(0.0)))

    gₚᴹ = (pₚᴹ*max(0, P - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    g_Dᴹ = (p_Dᴹ*max(0, D - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    gₚₒᴹ = (pₚₒᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    g_Zᴹ = (p_Zᴹ*max(0, Z - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    ∑gᴹ = gₚᴹ +  g_Dᴹ + gₚₒᴹ + g_Zᴹ #Sum grazing rates on each prey species for mesozooplankton
    
    return  ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ 
end

#GOC has variable sinking speed.
@inline function sinking_speed_of_GOC(z, zₑᵤ, zₘₓₗ, bgc)
    zₘₐₓ = max(abs(zₑᵤ), abs(zₘₓₗ)) 
    w_GOCᵐⁱⁿ = bgc.min_sinking_speed_of_GOC
    return w_GOCᵐⁱⁿ + (200/day - w_GOCᵐⁱⁿ)*(max(0, abs(z)-abs(zₘₐₓ)))/(5000) #41b
end

#Return flux feeding of mesozooplankton on POC and GOC, as well as a sum of flux feeding.
@inline function flux_feeding(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc) #eq29
    wₚₒ = bgc.sinking_speed_of_POC
    g_FF = bgc.flux_feeding_rate
    bₘ = bgc.temperature_sensitivity_term.M

    w_GOC = sinking_speed_of_GOC(z, zₑᵤ, zₘₓₗ, bgc)

    gₚₒ_FFᴹ = g_FF*(bₘ^T)*wₚₒ*POC #29a
    g_GOC_FFᴹ = g_FF*(bₘ^T)*w_GOC*GOC #29b
    ∑g_FFᴹ = g_GOC_FFᴹ + gₚₒ_FFᴹ
    return ∑g_FFᴹ, gₚₒ_FFᴹ, g_GOC_FFᴹ
end

#Gross growth efficiency is formulated to be called with either Z or M. However grazing on Z is only relevant for M, so pass zero when computing gross growth efficiency for Z.
@inline function nutrient_quality(gₚᴶ, g_Dᴶ, gₚₒᴶ, g_Zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    θᴺᶜ = bgc.NC_redfield_ratio
    θᶠᵉᶻ = bgc.FeC_ratio_of_zooplankton  #Assumed the same for both types of zooplankton

    ∑ᵢθᴺᴵgᵢᴶ = θᴺᶜ*gₚᴶ + θᴺᶜ*g_Dᴶ + θᴺᶜ*gₚₒᴶ + θᴺᶜ*g_Zᴹ
    ∑ᵢθᶠᵉᴵgᵢᴶ = nutrient_quota(Pᶠᵉ, P)*gₚᴶ + nutrient_quota(Dᶠᵉ, D)*g_Dᴶ + nutrient_quota(SFe, POC)*gₚₒᴶ + θᶠᵉᶻ*g_Zᴹ
    ∑ᵢgᵢᴶ = gₚᴶ + g_Dᴶ + gₚₒᴶ + g_Zᴹ
    
    return min(1, (∑ᵢθᴺᴵgᵢᴶ)/(θᴺᶜ*∑ᵢgᵢᴶ + eps(0.0)), (∑ᵢθᶠᵉᴵgᵢᴶ)/(θᶠᵉᶻ*∑ᵢgᵢᴶ + eps(0.0)))   #27a
end


@inline function growth_efficiency(eₘₐₓᴶ, σᴶ, gₚᴶ, g_Dᴶ, gₚₒᴶ, g_Zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    θᶠᵉᶻ = bgc.FeC_ratio_of_zooplankton  #Assumed the same for both types of zooplankton

    ∑ᵢθᶠᵉᴵgᵢᴶ = nutrient_quota(Pᶠᵉ, P)*gₚᴶ + nutrient_quota(Dᶠᵉ, D)*g_Dᴶ + nutrient_quota(SFe, POC)*gₚₒᴶ + θᶠᵉᶻ*g_Zᴹ
    ∑ᵢgᵢᴶ = gₚᴶ + g_Dᴶ + gₚₒᴶ + g_Zᴹ

    eₙᴶ = nutrient_quality(gₚᴶ, g_Dᴶ, gₚₒᴶ, g_Zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc) #27a

    return eₙᴶ*min(eₘₐₓᴶ, (1 - σᴶ)* (∑ᵢθᶠᵉᴵgᵢᴶ)/(θᶠᵉᶻ*∑ᵢgᵢᴶ + eps(0.0))) #27b
end


@inline function (bgc::PISCES)(::Val{:Z}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)    #args not correct
    #Parameters
    mᶻ = bgc.zooplankton_quadratic_mortality.Z
    b_Z = bgc.temperature_sensitivity_term.Z
    Kₘ = bgc.half_saturation_const_for_mortality
    rᶻ = bgc.zooplankton_linear_mortality.Z
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton.Z
    σᶻ = bgc.non_assimilated_fraction.Z

    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazing_Z(P, D, POC, T, bgc) 
    g_Zᴹ = grazing_M(P, D, Z, POC, T, bgc)[5]

    #Gross growth efficiency
    eᶻ = growth_efficiency(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    return (eᶻ*(gₚᶻ + g_Dᶻ + gₚₒᶻ)*Z - g_Zᴹ*M - mᶻ*(b_Z^T)*Z^2
        - rᶻ*(b_Z^T)*(concentration_limitation(Z, Kₘ) + 3*oxygen_conditions(O₂, bgc))*Z)   #24
end

@inline function (bgc::PISCES)(::Val{:M}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    mᴹ = bgc.zooplankton_quadratic_mortality.M
    bₘ = bgc.temperature_sensitivity_term.M
    rᴹ = bgc.zooplankton_linear_mortality.M
    Kₘ = bgc.half_saturation_const_for_mortality
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    σᴹ = bgc.non_assimilated_fraction.M

    #Grazing
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ  = grazing_M(P, D, Z, POC, T, bgc) 
    ∑g_FFᴹ = flux_feeding(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)[1]

    #Gross growth efficiency
    eᴹ =  growth_efficiency(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    return (eᴹ*(gₚᴹ + g_Dᴹ + gₚₒᴹ + ∑g_FFᴹ + g_Zᴹ)*M - mᴹ*(bₘ^T)*M^2 
        - rᴹ*(bₘ^T)*(concentration_limitation(M, Kₘ) + 3*oxygen_conditions(O₂, bgc))*M)   #28
end
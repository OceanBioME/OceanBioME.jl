@inline function get_grazingᶻ(P, D, POC, T, bgc) 
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

    grazing_arg = gₘₐₓᶻ*b_Z^T*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᶻ + pₚᶻ*P + p_Dᶻ*D + pₚₒᶻ*POC + eps(0.0)))

    gₚᶻ = (pₚᶻ*max(0, P - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    g_Dᶻ = (p_Dᶻ*max(0, D - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    gₚₒᶻ = (pₚₒᶻ*max(0, POC - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    
    ∑gᶻ= gₚᶻ + g_Dᶻ + gₚₒᶻ  #Sum grazing rates on each prey species for microzooplankton

    return ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ
end

@inline function get_grazingᴹ(P, D, Z, POC, T, bgc) 
    pₚᴹ = bgc.preference_for_nanophytoplankton.M
    p_Dᴹ = bgc.preference_for_diatoms.M
    pₚₒᴹ = bgc.preference_for_POC.M
    p_Zᴹ = bgc.preference_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    Fₜₕᵣₑₛₕᴹ = bgc.food_threshold_for_zooplankton.M
    gₘₐₓᴹ = bgc.max_grazing_rate.M
    K_Gᴹ = bgc.half_saturation_const_for_grazing.M
    bₘ = bgc.temperature_sensitivity_term.M
    
    F = pₚᴹ*max(0, P - Jₜₕᵣₑₛₕᴹ) + p_Dᴹ*max(0, D - Jₜₕᵣₑₛₕᴹ) + pₚₒᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ) + p_Zᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ)
    Fₗᵢₘ = max(0, F - min(0.5*F, Fₜₕᵣₑₛₕᴹ))

    grazing_arg =  gₘₐₓᴹ*bₘ^T*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᴹ + pₚᴹ*P + p_Dᴹ*D + pₚₒᴹ*POC + p_Zᴹ*Z + eps(0.0)))

    gₚᴹ = (pₚᴹ*max(0, P - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    g_Dᴹ = (p_Dᴹ*max(0, D - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    gₚₒᴹ = (pₚₒᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    g_Zᴹ = (p_Zᴹ*max(0, Z - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    ∑gᴹ = gₚᴹ +  g_Dᴹ + gₚₒᴹ + g_Zᴹ #Sum grazing rates on each prey species for mesozooplankton
    
    return  ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ 
end

@inline function get_w_GOC(z, zₑᵤ, zₘₓₗ, bgc)
    zₘₐₓ = max(zₑᵤ, zₘₓₗ) 
    w_GOCᵐⁱⁿ = bgc.min_sinking_speed_of_GOC
    return w_GOCᵐⁱⁿ + (200 - w_GOCᵐⁱⁿ)*(max(0, z-zₘₐₓ))/(5000) #41b
end

@inline function get_∑g_FFᴹ(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    wₚₒ = bgc.sinking_speed_of_POC
    g_FF = bgc.flux_feeding_rate
    bₘ = bgc.temperature_sensitivity_term.M

    w_GOC = get_w_GOC(z, zₑᵤ, zₘₓₗ, bgc)

    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC #29a
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC #29b
    return g_GOC_FFᴹ + gₚₒ_FFᴹ
end

# gross growth efficiency, defined for both but g_zᴹ and Z do not appear for eᶻ so have passed in as 0 
@inline function get_eₙᴶ(gₚᴶ, g_Dᴶ, gₚₒᴶ, g_Zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    θᴺᶜ = bgc.NC_redfield_ratio
    θᶠᵉᶜ = bgc.FeZ_redfield_ratio  #Assumed the same for both types of zooplankton

    ∑ᵢθᴺᴵgᵢᴶ = θᴺᶜ*gₚᴶ + θᴺᶜ*g_Dᴶ + θᴺᶜ*gₚₒᴶ + θᴺᶜ*g_Zᴹ
    ∑ᵢθᶠᵉᴵgᵢᴶ = θ(Pᶠᵉ, P)*gₚᴶ + θ(Dᶠᵉ, D)*g_Dᴶ + θ(SFe, POC)*gₚₒᴶ + θᶠᵉᶜ*g_Zᴹ
    ∑ᵢgᵢᴶ = gₚᴶ + g_Dᴶ + gₚₒᴶ + g_Zᴹ
    
    return min(1, (∑ᵢθᴺᴵgᵢᴶ)/(θᴺᶜ*∑ᵢgᵢᴶ), (∑ᵢθᶠᵉᴵgᵢᴶ)/(θᶠᵉᶜ*∑ᵢgᵢᴶ))   #27a
end


@inline function eᴶ(eₘₐₓᴶ, σᴶ, gₚᴶ, g_Dᴶ, gₚₒᴶ, g_Zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    θᶠᵉᶜ = bgc.FeZ_redfield_ratio  #Assumed the same for both types of zooplankton

    ∑ᵢθᶠᵉᴵgᵢᴶ = θ(Pᶠᵉ, P)*gₚᴶ + θ(Dᶠᵉ, D)*g_Dᴶ + θ(SFe, POC)*gₚₒᴶ + θᶠᵉᶜ*g_Zᴹ
    ∑ᵢgᵢᴶ = gₚᴶ + g_Dᴶ + gₚₒᴶ + g_Zᴹ

    eₙᴶ = get_eₙᴶ(gₚᴶ, g_Dᴶ, gₚₒᴶ, g_Zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc) #27a

    return eₙᴶ*min(eₘₐₓᴶ, (1 - σᴶ)* (∑ᵢθᶠᵉᴵgᵢᴶ)/(θᶠᵉᶜ*∑ᵢgᵢᴶ)) #27b
end


@inline function (bgc::PISCES)(::Val{:Z}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust)    #args not correct
    mᶻ = bgc.zooplankton_quadratic_mortality.Z
    b_Z = bgc.temperature_sensitivity_term.Z
    Kₘ = bgc.half_saturation_const_for_mortality
    rᶻ = bgc.zooplankton_linear_mortality.Z
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton.Z
    σᶻ = bgc.non_assimilated_fraction.Z

    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = get_grazingᶻ(P, D, POC, T, bgc) 
    g_Zᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)[5]

    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    return eᶻ*(gₚᶻ + g_Dᶻ + gₚₒᶻ)*Z - g_Zᴹ*M - mᶻ*b_Z^T*Z^2 - rᶻ*b_Z^T*(K_mondo(Z, Kₘ) + 3*ΔO₂(O₂, bgc))*Z   #24
end

@inline function (bgc::PISCES)(::Val{:M}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust)
    mᴹ = bgc.zooplankton_quadratic_mortality.M
    bₘ = bgc.temperature_sensitivity_term.M
    rᴹ = bgc.zooplankton_linear_mortality.M
    Kₘ = bgc.half_saturation_const_for_mortality

    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    σᴹ = bgc.non_assimilated_fraction.M

    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ  = get_grazingᴹ(P, D, Z, POC, T, bgc) 

    ∑g_FFᴹ = get_∑g_FFᴹ(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    return eᴹ*(gₚᴹ + g_Dᴹ + gₚₒᴹ + ∑g_FFᴹ + g_Zᴹ)*M - mᴹ*bₘ^T*M^2 - rᴹ*bₘ^T*(K_mondo(M, Kₘ) + 3*ΔO₂(O₂, bgc))*M   #28
end
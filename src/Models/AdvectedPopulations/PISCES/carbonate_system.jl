#This document contains functions for:
    #Forcing for DIC.
    #Forcing for Alk.

@inline function (bgc::PISCES)(::Val{:DIC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³)
    #Parameters
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    eₘₐₓᶻ = bgc. max_growth_efficiency_of_zooplankton.Z
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    αᴾ= bgc.initial_slope_of_PI_curve.P
    αᴰ = bgc.initial_slope_of_PI_curve.D

    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)

    bFe = Fe

    zₘₐₓ = max(zₘₓₗ, zₑᵤ)
    
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = get_grazingᶻ(P, D, POC, T, bgc)
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)
    ∑g_FFᴹ = get_∑g_FFᴹ(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)

    #PAR
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = get_PARᴾ(PAR¹, PAR², PAR³, bgc)
    PARᴰ = get_PARᴰ(PAR¹, PAR², PAR³, bgc)

    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    Bact = get_Bact(zₘₐₓ, z, Z, M)

    #Growth rates for phytoplankton
    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc)
    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)

    return γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M + γᴹ*Rᵤₚ(M, T, bgc) + get_Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc) + get_Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc) + λ_CaCO₃¹(CaCO₃, bgc, Ω)*CaCO₃ - P_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, D, Z, M, POC, T, PAR, zₘₓₗ, z, bgc) - μᴰ*D - μᴾ*P #eq59
end

@inline function (bgc::PISCES)(::Val{:Alk}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³) # eq59
    #Parameters
    θᴺᶜ = bgc.NC_redfield_ratio
    rₙₒ₃¹ = bgc. CN_ratio_of_denitrification
    rₙₕ₄¹ = bgc.CN_ratio_of_ammonification
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    λₙₕ₄ = bgc.max_nitrification_rate
    eₘₐₓᶻ = bgc. max_growth_efficiency_of_zooplankton.Z
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M

    bFe = Fe

    zₘₐₓ = max(zₘₓₗ, zₑᵤ)

    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)

    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = get_PARᴾ(PAR¹, PAR², PAR³, bgc)
    PARᴰ = get_PARᴰ(PAR¹, PAR², PAR³, bgc)

    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = get_grazingᶻ(P, D, POC, T, bgc)
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)
    ∑g_FFᴹ = get_∑g_FFᴹ(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)

    Bact = get_Bact(zₘₐₓ, z, Z, M)

    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc) 
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
   
    #Uptake rates of nitrogen and ammonium
    μₙₒ₃ᴾ = get_μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₒ₃ᴰ = get_μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)
    μₙₕ₄ᴾ = get_μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₕ₄ᴰ = get_μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)

    return θᴺᶜ*get_Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc) + θᴺᶜ*(rₙₒ₃¹ + 1)*get_Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc) + θᴺᶜ*γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + θᴺᶜ*γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ + θᴺᶜ*γᴹ*Rᵤₚ(M, T, bgc))*M + θᴺᶜ*μₙₒ₃ᴾ*P + θᴺᶜ*μₙₒ₃ᴰ*D + θᴺᶜ*N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR, bgc) + 2*λ_CaCO₃¹(CaCO₃, bgc, Ω)*CaCO₃ + θᴺᶜ*ΔO₂(O₂, bgc)*(rₙₕ₄¹ - 1)*λₙₕ₄*NH₄ - θᴺᶜ*μₙₕ₄ᴾ*P - θᴺᶜ*μₙₕ₄ᴰ*D- 2*θᴺᶜ*Nitrif(NH₄, O₂, λₙₕ₄, PAR, bgc) - 2*P_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, D, Z, M, POC, T, PAR, zₘₓₗ, z, bgc)
end
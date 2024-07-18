#This document contains functions for:
    #Forcing for DIC.
    #Forcing for Alk.

@inline function (bgc::PISCES)(::Val{:DIC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust)
    #Parameters
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    eₘₐₓᶻ = bgc. max_growth_efficiency_of_zooplankton.Z
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M

    bFe = Fe
    
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazingᶻ(P, D, POC, T, bgc)
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ = grazingᴹ(P, D, Z, POC, T, bgc)
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    
    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z, bgc)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    #Growth rates for phytoplankton
    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc)
    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)

    return γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M + γᴹ*Rᵤₚᴹ(M, T) 
    + Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc) + Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc) 
    + λ_CaCO₃¹(CaCO₃, bgc)*CaCO₃ - P_CaCO₃(P, Z, M, T, PAR, zₘₓₗ, z, bgc) - μᴰ*D - μᴾ*P #eq59
end

@inline function (bgc::PISCES)(::Val{:Alk}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust) # eq59
    #Parameters
    θᴺᶜ = bgc.NC_redfield_ratio
    rₙₒ₃¹ = bgc. CN_ratio_of_denitrification
    rₙₕ₄¹ = bgc.CN_ratio_of_ammonification
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    λₙₕ₄ = bgc.max_nitrification_rate

    bFe = Fe

    #Grazing
    grazingᶻ = grazingᶻ(P, D, POC, T, bgc)
    grazingᴹ = grazingᴹ(P, D, Z, POC, T, bgc)
    ∑gᶻ = grazingᶻ[1]
    ∑gᴹ = grazingᴹ[1]
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)

    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z, bgc)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
   
    #Uptake rates of nitrogen and ammonium
    μₙₒ₃ᴾ = μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, bgc)
    μₙₒ₃ᴰ = μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, bgc)
    μₙₕ₄ᴾ = μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ)
    μₙₕ₄ᴰ = μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ)

    return θᴺᶜ*Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact) + θᴺᶜ*(rₙₒ₃¹ + 1)*Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact) 
    + θᴺᶜ*γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + θᴺᶜ*γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ + θᴺᶜ*γᴹ*Rᵤₚ(M, T))*M + θᴺᶜ*μₙₒ₃ᴾ*P + θᴺᶜ*μₙₒ₃ᴰ*D 
    + θᴺᶜ*N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR) + 2*λ_CaCO₃¹(CaCO₃)*CaCO₃ + θᴺᶜ*ΔO₂(O₂)*(rₙₕ₄¹ - 1)*λₙₕ₄*NH₄ 
    - θᴺᶜ*μₙₕ₄ᴾ*P - θᴺᶜ*μₙₕ₄ᴰ*D- 2*θᴺᶜ*Nitrif(NH₄, O₂, λₙₕ₄, PAR) - 2*P_CaCO₃(P, Z, M, T, PAR, zₘₓₗ)
end
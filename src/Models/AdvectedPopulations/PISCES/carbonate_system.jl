#This document contains functions for:
    #Forcing for DIC.
    #Forcing for Alk.

@inline function (pisces::PISCES)(::Val{:DIC}, x, y, z, t, P, D, Z, M, CaCO₃, T, PAR, zₘₓₗ)
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    
    #Grazing
    grazingᶻ = grazingᶻ(P, D, POC, T)
    grazingᴹ = grazingᴹ(P, D, Z, POC, T)
    ∑gᶻ = grazingᶻ[1]
    ∑gᴹ = grazingᴹ[1]
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)

    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC)

    #Growth rates for phytoplankton
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ)
    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)

    return γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M + γᴹ*Rᵤₚᴹ(M, T) + Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact) + Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact) + λ_CaCO₃¹()*CaCO₃ - P_CaCO₃(P, Z, M, T, PAR, zₘₓₗ) - μᴰ*D - μᴾ*P #eq59
end

@inline function (pisces::PISCES)(::Val{:Alk}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅) # eq59
    
    θᴺᶜ = bgc.NC_redfield_ratio
    rₙₒ₃¹ = bgc. CN_ratio_of_denitrification
    rₙₕ₄¹ = bgc.CN_ratio_of_ammonification
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    λₙₕ₄ = bgc.max_nitrification_rate

    #Grazing
    grazingᶻ = grazingᶻ(P, D, POC, T)
    grazingᴹ = grazingᴹ(P, D, Z, POC, T)
    ∑gᶻ = grazingᶻ[1]
    ∑gᴹ = grazingᴹ[1]
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)

    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC)
   
    #Uptake rates of nitrogen and ammonium
    μₙₒ₃ᴾ = μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ)
    μₙₒ₃ᴰ = μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ)
    μₙₕ₄ᴾ = μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ)
    μₙₕ₄ᴰ = μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ)

    return θᴺᶜ*Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact) + θᴺᶜ*(rₙₒ₃¹ + 1)*Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact) + θᴺᶜ*γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + θᴺᶜ*γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ + θᴺᶜ*γᴹ*Rᵤₚ(M, T))*M + θᴺᶜ*μₙₒ₃ᴾ*P + θᴺᶜ*μₙₒ₃ᴰ*D + θᴺᶜ*N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR) + 2*λ_CaCO₃¹()*CaCO₃ + θᴺᶜ*ΔO₂(O₂)*(rₙₕ₄¹ - 1)*λₙₕ₄*NH₄ - θᴺᶜ*μₙₕ₄ᴾ*P - θᴺᶜ*μₙₕ₄ᴰ*D- 2*θᴺᶜ*Nitrif(NH₄, O₂, λₙₕ₄, PAR) - 2*P_CaCO₃(P, Z, M, T, PAR, zₘₓₗ)
end
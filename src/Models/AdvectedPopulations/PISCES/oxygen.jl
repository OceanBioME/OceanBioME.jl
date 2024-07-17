#This document contains functions for:
    #O₂ forcing (eq83)

@inline function (pisces::PISCES)(::Val{:O₂}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust) 

    O₂ᵘᵗ = bgc.OC_for_ammonium_based_processes
    O₂ⁿⁱᵗ = bgc.OC_ratio_of_nitrification
    γᶻ = bgc.excretion_as_DOM.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᶻ = bgc.non_assimilated_fraction.Z
    σᴹ = bgc.non_assimilated_fraction.M

    #L_day
    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)
    
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
    
    return O₂ᵘᵗ*(μₙₕ₄ᴾ*P + μₙₕ₄ᴰ*D) + (O₂ᵘᵗ + O₂ⁿⁱᵗ)*(μₙₒ₃ᴾ*P + μₙₒ₃ᴰ*D) + O₂ⁿⁱᵗ*N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR) - O₂ᵘᵗ*γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z - O₂ᵘᵗ*γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M - O₂ᵘᵗ*γᴹ*Rᵤₚᴹ(M, T) - O₂ᵘᵗ*Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact) - O₂ⁿⁱᵗ*Nitrif(NH₄, O₂, λₙₕ₄, PAR)
end
#This document contains functions for:
    #O₂ forcing (eq83)

@inline function (bgc::PISCES)(::Val{:O₂}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust) 

    O₂ᵘᵗ = bgc.OC_for_ammonium_based_processes
    O₂ⁿⁱᵗ = bgc.OC_ratio_of_nitrification
    γᶻ = bgc.excretion_as_DOM.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᶻ = bgc.non_assimilated_fraction.Z
    σᴹ = bgc.non_assimilated_fraction.M
    λₙₕ₄ = bgc.max_nitrification_rate
    bFe = Fe

    #L_day
    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)

    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = get_PARᴾ(PAR¹, PAR², PAR³, bgc)
    PARᴰ = get_PARᴰ(PAR¹, PAR², PAR³, bgc)
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton.Z
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = get_grazingᶻ(P, D, POC, T, bgc) 
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ  = get_grazingᴹ(P, D, Z, POC, T, bgc) 
  
    ∑g_FFᴹ = get_∑g_FFᴹ(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    #g_Z not called
    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    zₘₐₓ = max(zₑᵤ, zₘₓₗ) #35a
    Bact = get_Bact(zₘₐₓ, z, Z, M)
   
    #Uptake rates of nitrogen and ammonium
    μₙₒ₃ᴾ = get_μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₒ₃ᴰ = get_μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)
    μₙₕ₄ᴾ = get_μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₕ₄ᴰ = get_μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)
    
    return O₂ᵘᵗ*(μₙₕ₄ᴾ*P + μₙₕ₄ᴰ*D) + (O₂ᵘᵗ + O₂ⁿⁱᵗ)*(μₙₒ₃ᴾ*P + μₙₒ₃ᴰ*D) + O₂ⁿⁱᵗ*N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR, bgc) - O₂ᵘᵗ*γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z - O₂ᵘᵗ*γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M - O₂ᵘᵗ*γᴹ*Rᵤₚ(M, T, bgc) - O₂ᵘᵗ*get_Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc) - O₂ⁿⁱᵗ*Nitrif(NH₄, O₂, λₙₕ₄, PAR, bgc)
end
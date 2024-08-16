#This document contains functions for:
    #O₂ forcing (eq83)

@inline function (bgc::PISCES)(::Val{:O₂}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³) 

    O₂ᵘᵗ = bgc.OC_for_ammonium_based_processes
    O₂ⁿⁱᵗ = bgc.OC_ratio_of_nitrification
    γᶻ = bgc.excretion_as_DOM.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᶻ = bgc.non_assimilated_fraction.Z
    σᴹ = bgc.non_assimilated_fraction.M
    λₙₕ₄ = bgc.max_nitrification_rate
    bFe = Fe
    θᴺᶜ = bgc.NC_redfield_ratio
    #L_day
    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = latitude(ϕ₀, y)
    L_day = day_length(ϕ, t, L_day_param)

    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = P_PAR(PAR¹, PAR², PAR³, bgc)
    PARᴰ = D_PAR(PAR¹, PAR², PAR³, bgc)
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton.Z
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazing_Z(P, D, POC, T, bgc) 
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ  = grazing_M(P, D, Z, POC, T, bgc) 
    ∑g_FFᴹ = flux_feeding(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)[1]
    #g_Z not called
    
    #Gross growth efficiency
    eᶻ = growth_efficiency(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eᴹ =  growth_efficiency(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    zₘₐₓ = max(abs(zₑᵤ), abs(zₘₓₗ)) #35a
    Bact = bacterial_biomass(zₘₐₓ, z, Z, M)
   
    #Uptake rates of nitrogen and ammonium
    μₙₒ₃ᴾ = uptake_rate_nitrate_P(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₒ₃ᴰ = uptake_rate_nitrate_D(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)
    μₙₕ₄ᴾ = uptake_rate_ammonium_P(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₕ₄ᴰ = uptake_rate_ammonium_D(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)
    
    return (O₂ᵘᵗ*(μₙₕ₄ᴾ*P + μₙₕ₄ᴰ*D) + (O₂ᵘᵗ + O₂ⁿⁱᵗ)*(μₙₒ₃ᴾ*P + μₙₒ₃ᴰ*D) + O₂ⁿⁱᵗ*(1/θᴺᶜ)*N_fixation(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR, bgc) 
            - O₂ᵘᵗ*γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z - O₂ᵘᵗ*γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M - O₂ᵘᵗ*γᴹ*upper_respiration(M, T, bgc)
            - O₂ᵘᵗ*oxic_remineralization(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc) - O₂ⁿⁱᵗ*(1/θᴺᶜ)*nitrification(NH₄, O₂, λₙₕ₄, PAR, bgc))
end
#This document contains functions for:
    #PO₄ forcing (eq59), multiplied by redfield ratio to return in μmolP/L

@inline function (bgc::PISCES)(::Val{:PO₄}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃) 
    #Parameters
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    αᴾ = bgc.initial_slope_of_PI_curve.P
    αᴰ = bgc.initial_slope_of_PI_curve.D
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton.Z
    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton.M
    θᴾᶜ = bgc.PC_redfield_ratio

    bFe = Fe
    
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazing_Z(P, D, POC, T, bgc) 
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ  = grazing_M(P, D, Z, POC, T, bgc) 
    ∑g_FFᴹ = flux_feeding(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)[1]

    #Gross growth efficiency
    eᶻ = growth_efficiency(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eᴹ =  growth_efficiency(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    #Bacteria
    zₘₐₓ = max(abs(zₑᵤ), abs(zₘₓₗ)) #35a
    Bact = bacterial_biomass(zₘₐₓ, z, Z, M)

    #Growth rates for phytoplankton
    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = latitude(ϕ₀, y)
    L_day = day_length(ϕ, t, L_day_param)
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = P_PAR(PAR₁, PAR₂, PAR₃, bgc)
    PARᴰ = D_PAR(PAR₁, PAR₂, PAR₃, bgc)

    Lₗᵢₘᴾ = P_nutrient_limitation(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    μᴾ = phytoplankton_growth_rate(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc)
    μᴰ = phytoplankton_growth_rate(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)

    return (θᴾᶜ*(γᶻ*(1-eᶻ-σᶻ)*∑gᶻ*Z + γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M + γᴹ*upper_respiration(M, T, bgc) 
            + oxic_remineralization(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc) + denitrification(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc) 
            - μᴾ*P  - μᴰ*D)) #eq59
end
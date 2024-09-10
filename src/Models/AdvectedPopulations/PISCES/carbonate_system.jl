#This document contains functions for:
    #Forcing for DIC.
    #Forcing for Alk.

#DIC is significant as required by phytoplankton for photosynthesis, and by calcifying organisms for calcite shells.

@inline function (bgc::PISCES)(::Val{:DIC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR², PAR³)
    #Parameters
    γᶻ, γᴹ = bgc.excretion_as_DOM
    σᶻ, σᴹ = bgc.non_assimilated_fraction
    eₘₐₓᶻ, eₘₐₓᴹ = bgc. max_growth_efficiency_of_zooplankton
    αᴾ, αᴰ = bgc.initial_slope_of_PI_curve

    bFe = Fe
    
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazing_Z(P, D, POC, T, bgc)
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ = grazing_M(P, D, Z, POC, T, bgc)
    ∑g_FFᴹ = flux_feeding(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)[1]

    #Gross growth efficiency
    eᶻ = growth_efficiency(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eᴹ =  growth_efficiency(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    #Growth rates for phytoplankton
    φ = bgc.latitude
    φ = latitude(φ, y)


    L_day = day_length(ϕ, t)
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = P_PAR(PAR₁, PAR², PAR³, bgc)
    PARᴰ = D_PAR(PAR₁, PAR², PAR³, bgc)

    Lₗᵢₘᴾ = P_nutrient_limitation(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[1]
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    μᴾ = phytoplankton_growth_rate(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ, bgc)
    μᴰ = phytoplankton_growth_rate(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)

    #Bacteria
    zₘₐₓ = max(abs(zₘₓₗ), abs(zₑᵤ))
    Bact = bacterial_biomass(zₘₐₓ, z, Z, M)

    return (γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M + 
    γᴹ*upper_respiration(M, T, bgc) + 
    oxic_remineralization(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc) + 
    denitrification(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc) + 
    dissolution_of_calcite(CaCO₃, bgc, Ω)*CaCO₃ - 
    production_of_sinking_calcite(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, D, Z, M, POC, T, PAR, zₘₓₗ, z, bgc) - 
    μᴰ*D - μᴾ*P) #eq59
end

@inline function (bgc::PISCES)(::Val{:Alk}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR², PAR³) # eq59
    #Parameters
    θᴺᶜ = bgc.NC_redfield_ratio
    rₙₒ₃¹ = bgc. CN_ratio_of_denitrification
    rₙₕ₄¹ = bgc.CN_ratio_of_ammonification
    γᶻ, γᴹ = bgc.excretion_as_DOM
    σᶻ, σᴹ = bgc.non_assimilated_fraction
    eₘₐₓᶻ, eₘₐₓᴹ = bgc. max_growth_efficiency_of_zooplankton
    λₙₕ₄ = bgc.max_nitrification_rate

    bFe = Fe

    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazing_Z(P, D, POC, T, bgc)
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ = grazing_M(P, D, Z, POC, T, bgc)
    ∑g_FFᴹ = flux_feeding(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)[1]

    #Gross growth efficiency
    eᶻ = growth_efficiency(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc) 
    eᴹ =  growth_efficiency(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    #Uptake rates of nitrogen and ammonium
    φ = bgc.latitude
    φ = latitude(φ, y)


    L_day = day_length(ϕ, t)
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    PARᴾ = P_PAR(PAR₁, PAR², PAR³, bgc)
    PARᴰ = D_PAR(PAR₁, PAR², PAR³, bgc)

    μₙₒ₃ᴾ = uptake_rate_nitrate_P(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₒ₃ᴰ = uptake_rate_nitrate_D(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)
    μₙₕ₄ᴾ = uptake_rate_ammonium_P(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ, Si̅, bgc)
    μₙₕ₄ᴰ = uptake_rate_ammonium_D(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ, Si̅, bgc)

    #Bacteria
    zₘₐₓ = max(abs(zₘₓₗ), abs(zₑᵤ))
    Bact = bacterial_biomass(zₘₐₓ, z, Z, M)

    return (θᴺᶜ*oxic_remineralization(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact, bgc)
          + θᴺᶜ*(rₙₒ₃¹ + 1)*denitrification(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact, bgc)
          + θᴺᶜ*γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + θᴺᶜ*γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ 
          + θᴺᶜ*γᴹ*upper_respiration(M, T, bgc))*M + θᴺᶜ*μₙₒ₃ᴾ*P + θᴺᶜ*μₙₒ₃ᴰ*D 
          + N_fixation(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR, bgc)
          + 2*dissolution_of_calcite(CaCO₃, bgc, Ω)*CaCO₃ + θᴺᶜ*oxygen_conditions(O₂, bgc)*(rₙₕ₄¹ - 1)*λₙₕ₄*NH₄
          - θᴺᶜ*μₙₕ₄ᴾ*P - θᴺᶜ*μₙₕ₄ᴰ*D- 2*θᴺᶜ*nitrification(NH₄, O₂, λₙₕ₄, PAR, bgc) 
          - 2*production_of_sinking_calcite(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, D, Z, M, POC, T, PAR, zₘₓₗ, z, bgc)) #eq81
end

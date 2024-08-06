
# Using simple chemistry model. 
# This document contains functions for the following:
    # Fe¹ (eq65), dissolved free inorganic iron
    # Cgfe1, Cgfe2, Aggfe, Bactfe (eqs 61, 62, 63)
    # Forcing for Fe (eq60)

@inline function get_Fe¹(Fe, DOC, T)
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6) # bgc.total_concentration_of_iron_ligands
    K_eqᶠᵉ = 10^(16.27 - 1565.7/max(T + 273.15, 5)) #check this value
    Δ = 1 +  K_eqᶠᵉ*Lₜ -  K_eqᶠᵉ*Fe

    return (-Δ + sqrt(Δ^2 + 4*K_eqᶠᵉ*Fe))/(2*K_eqᶠᵉ + eps(0.0)) #eq65
end

@inline function Cgfe1(sh, Fe, POC, DOC, T, bgc)
    a₁ = bgc.aggregation_rate_of_DOC_to_POC_1
    a₂ = bgc.aggregation_rate_of_DOC_to_POC_2
    a₄ = bgc.aggregation_rate_of_DOC_to_POC_4
    a₅ = bgc.aggregation_rate_of_DOC_to_POC_5
   
    FeL = Fe - get_Fe¹(Fe, DOC, T) #eq64
    Fe_coll = 0.5*FeL
    return ((a₁*DOC + a₂*POC)*sh+a₄*POC + a₅*DOC)*Fe_coll
end

@inline function Cgfe2(sh, Fe, T, DOC, GOC, bgc)
    a₃ = bgc.aggregation_rate_of_DOC_to_GOC_3
    FeL = Fe - get_Fe¹(Fe, DOC, T)
    Fe_coll = 0.5*FeL
    return a₃*GOC*sh*Fe_coll
end

@inline function Aggfe(Fe, DOC, T, bgc)
    λᶠᵉ = 1e-3 * bgc.slope_of_scavenging_rate_of_iron #parameter not defined in parameter list. Assumed scaled version λ_Fe to fit dimensions of Fe¹.
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6)
    return 1000*λᶠᵉ*max(0, Fe - Lₜ)*get_Fe¹(Fe, DOC, T)
end

@inline function get_Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc)
    K_Feᴮ¹ = bgc.Fe_half_saturation_const_for_PLACEHOLDER
    θₘₐₓᶠᵉᵇᵃᶜᵗ = bgc.max_FeC_ratio_of_bacteria
    Bact = get_Bact(zₘₐₓ, z, Z, M) 
    Lₗᵢₘᵇᵃᶜᵗ = Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe, bgc)[2]
    return μₘₐₓ⁰*fₚ(T, bgc)*Lₗᵢₘᵇᵃᶜᵗ*θₘₐₓᶠᵉᵇᵃᶜᵗ*Fe*Bact/(K_Feᴮ¹ + Fe + eps(0.0)) 
end

@inline function (bgc::PISCES)(::Val{:Fe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³) #eq60
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    δᴾ = bgc.exudation_of_DOC.P
    δᴰ = bgc.exudation_of_DOC.D
    θᶠᵉᶻ = bgc.FeC_ratio_of_zooplankton
    μₘₐₓ⁰ = bgc.growth_rate_at_zero
    θₘₐₓᶠᵉᴾ = bgc.max_iron_quota.P
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton.P
    K_Feᴾᶠᵉᵐⁱⁿ = bgc.min_half_saturation_const_for_iron_uptake.P
    Pₘₐₓ = bgc.threshold_concentration_for_size_dependency.P
    θₘₐₓᶠᵉᴰ = bgc.max_iron_quota.D
    Sᵣₐₜᴰ = bgc.size_ratio_of_phytoplankton.D
    K_Feᴰᶠᵉᵐⁱⁿ = bgc.min_half_saturation_const_for_iron_uptake.D
    Dₘₐₓ = bgc.threshold_concentration_for_size_dependency.D
    g_FF = bgc.flux_feeding_rate
    bₘ = bgc.temperature_sensitivity_term.M
    wₚₒ = bgc.sinking_speed_of_POC
    w_GOCᵐⁱⁿ = bgc.min_sinking_speed_of_GOC

    bFe = Fe

    L_Feᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[6]
    L_Feᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[6]
    
    sh = get_sh(z, zₘₓₗ)
    
    λₚₒ¹ = λ¹(T, O₂, bgc)

    μᴾᶠᵉ = μᴵᶠᵉ(P, Pᶠᵉ, θₘₐₓᶠᵉᴾ, Sᵣₐₜᴾ, K_Feᴾᶠᵉᵐⁱⁿ, Pₘₐₓ, L_Feᴾ, bFe, T, bgc)
    μᴰᶠᵉ = μᴵᶠᵉ(D, Dᶠᵉ, θₘₐₓᶠᵉᴰ, Sᵣₐₜᴰ, K_Feᴰᶠᵉᵐⁱⁿ, Dₘₐₓ, L_Feᴰ, bFe, T, bgc)
    zₘₐₓ = max(zₑᵤ, zₘₓₗ)
    #Iron quotas
    θᶠᵉᴾ = θ(Pᶠᵉ, P)
    θᶠᵉᴰ = θ(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = θ(SFe, POC)
    θᶠᵉᴳᴼᶜ = θ(BFe, GOC)
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = get_grazingᶻ(P, D, POC, T, bgc)
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)
    ∑g_FFᴹ = get_∑g_FFᴹ(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    
    ∑θᶠᵉⁱgᵢᶻ = θᶠᵉᴾ*gₚᶻ + θᶠᵉᴰ*g_Dᶻ + θᶠᵉᴾᴼᶜ*gₚₒᶻ #over P, D, POC
    ∑θᶠᵉⁱgᵢᴹ = θᶠᵉᴾ*gₚᴹ + θᶠᵉᴰ*g_Dᴹ + θᶠᵉᴾᴼᶜ*gₚₒᴹ + θᶠᵉᶻ*g_Zᴹ #graze on P, D, POC, Z 

    Bactfe = get_Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc)

    gₚₒ_FF = g_FF*bₘ^T*wₚₒ*POC#
    w_GOC = w_GOCᵐⁱⁿ + (200 - w_GOCᵐⁱⁿ)*(max(0, z-zₘₐₓ))/(5000) #41b
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC 

    #Gross growth efficiency
    eₙᶻ = get_eₙᴶ(gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eₙᴹ = get_eₙᴶ(gₚᴹ, g_Dᴹ, gₚₒᴹ, g_Zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)

    #println("Sum of positive terms in Iron is ",  max(0, (1-σᶻ)*(∑θᶠᵉⁱgᵢᶻ/(∑gᶻ + eps(0.0)) - eₙᶻ*θᶠᵉᶻ))*∑gᶻ*Z + max(0, (1-σᴹ)*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FF + θᶠᵉᴳᴼᶜ*g_GOC_FFᴹ )/(∑gᴹ+∑g_FFᴹ + eps(0.0)) - eₙᴹ*θᶠᵉᶻ)*(∑gᴹ+∑g_FFᴹ)*M + γᴹ*θᶠᵉᶻ*Rᵤₚ(M, T, bgc) + λₚₒ¹*SFe)
    #println("term a1 = $(max(0, (1-σᶻ)*(∑θᶠᵉⁱgᵢᶻ/(∑gᶻ + eps(0.0)) - eₙᶻ*θᶠᵉᶻ))*∑gᶻ*Z), term b = $(max(0, (1-σᴹ)*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FF + θᶠᵉᴳᴼᶜ*g_GOC_FFᴹ )/(∑gᴹ+∑g_FFᴹ + eps(0.0)) - eₙᴹ*θᶠᵉᶻ)*(∑gᴹ+∑g_FFᴹ)*M ), term c = $( γᴹ*θᶠᵉᶻ*Rᵤₚ(M, T, bgc)), term d = $(λₚₒ¹*SFe)")
    #println("Sum of negative terms is ", (1 - δᴾ)*μᴾᶠᵉ*P + (1 - δᴰ)*μᴰᶠᵉ*D + Scav(POC, GOC, CaCO₃, PSi, D_dust, DOC, T, Fe, bgc) + Cgfe1(sh, Fe, POC, DOC, T, bgc) + Cgfe2(sh, Fe, T, DOC, GOC, bgc) + Aggfe(Fe, DOC, T, bgc) + Bactfe)
    #println("term 1 = $((1 - δᴾ)*μᴾᶠᵉ*P ), term 2 = $((1 - δᴰ)*μᴰᶠᵉ*D), term 3 = $(Scav(POC, GOC, CaCO₃, PSi, D_dust, DOC, T, Fe, bgc)), term 4 = $(Cgfe1(sh, Fe, POC, DOC, T, bgc)), term 5 = $(Cgfe2(sh, Fe, T, DOC, GOC, bgc)), term 6 = $(Aggfe(Fe, DOC, T, bgc)), term 6 = $(Bactfe)")
    #println("--")
    #println("Total change is ",max(0, (1-σᶻ)*(∑θᶠᵉⁱgᵢᶻ/(∑gᶻ + eps(0.0)) - eₙᶻ*θᶠᵉᶻ))*∑gᶻ*Z + max(0, (1-σᴹ)*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FF + θᶠᵉᴳᴼᶜ*g_GOC_FFᴹ )/(∑gᴹ+∑g_FFᴹ + eps(0.0)) - eₙᴹ*θᶠᵉᶻ)*(∑gᴹ+∑g_FFᴹ)*M + γᴹ*θᶠᵉᶻ*Rᵤₚ(M, T, bgc) + λₚₒ¹*SFe - (1 - δᴾ)*μᴾᶠᵉ*P - (1 - δᴰ)*μᴰᶠᵉ*D - Scav(POC, GOC, CaCO₃, PSi, D_dust, DOC, T, Fe, bgc) - Cgfe1(sh, Fe, POC, DOC, T, bgc) - Cgfe2(sh, Fe, T, DOC, GOC, bgc) - Aggfe(Fe, DOC, T, bgc) - Bactfe)
    #println("-------------------------------------")

    return max(0, (1-σᶻ)*(∑θᶠᵉⁱgᵢᶻ/(∑gᶻ + eps(0.0)) - eₙᶻ*θᶠᵉᶻ))*∑gᶻ*Z + max(0, (1-σᴹ)*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FF + θᶠᵉᴳᴼᶜ*g_GOC_FFᴹ )/(∑gᴹ+∑g_FFᴹ + eps(0.0)) - eₙᴹ*θᶠᵉᶻ)*(∑gᴹ+∑g_FFᴹ)*M + γᴹ*θᶠᵉᶻ*Rᵤₚ(M, T, bgc) + λₚₒ¹*SFe - (1 - δᴾ)*μᴾᶠᵉ*P - (1 - δᴰ)*μᴰᶠᵉ*D - Scav(POC, GOC, CaCO₃, PSi, D_dust, DOC, T, Fe, bgc) - Cgfe1(sh, Fe, POC, DOC, T, bgc) - Cgfe2(sh, Fe, T, DOC, GOC, bgc) - Aggfe(Fe, DOC, T, bgc) - Bactfe
end

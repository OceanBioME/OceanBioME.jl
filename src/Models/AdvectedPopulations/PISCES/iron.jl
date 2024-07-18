
# Using simple chemistry model. 
# This document contains functions for the following:
    # Fe¹ (eq65), dissolved free inorganic iron
    # Cgfe1, Cgfe2, Aggfe, Bactfe (eqs 61, 62, 63)
    # Forcing for Fe (eq60)

@inline function Fe¹(Fe, DOC, T)
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6) # bgc.total_concentration_of_iron_ligands
    K_eqᶠᵉ = 10^(16.27 - 1565.7/max(T, 5)) #check this value
    Δ = 1 +  K_eqᶠᵉ(T)*Lₜ -  K_eqᶠᵉ(T)*Fₑ

    return (-Δ + sqrt(Δ^2 + 4*K_eqᶠᵉ(T)*Fe))/(2*K_eqᶠᵉ(T)) #eq65
end

@inline function Cgfe1(sh, Fe, POC, DOC, T, bgc)
    a₁ = bgc.aggregation_rate_of_DOC_to_POC_1
    a₂ = bgc.aggregation_rate_of_DOC_to_POC_2
    a₄ = bgc.aggregation_rate_of_DOC_to_POC_4
    a₅ = bgc.aggregation_rate_of_DOC_to_POC_5
   
    FeL = Fe - Fe¹(DOC, T, Fe) #eq64
    Fe_coll = 0.5*FeL
    return ((a₁*DOC + a₂*POC)*sh+a₄*POC + a₅*DOC)*Fe_coll
end

@inline function Cgfe2(sh, Fe, T, DOC, GOC, bgc)
    a₃ = bgc.aggregation_rate_of_DOC_to_GOC_3
    FeL = Fe - Fe¹(Fe, DOC, T)
    Fe_coll = 0.5*FeL
    return a₃*GOC*sh*Fe_coll
end

@inline function Aggfe(Fe, DOC, T, bgc)
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6)
    return 1000*λ_Fe*max(0, Fe - Lₜ)*Fe¹(DOC, T, Fe)
end

@inline function Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc)
    K_Feᴮ¹ = bgc.Fe_half_saturation_const_for_PLACEHOLDER
    θₘₐₓᶠᵉᵇᵃᶜᵗ = bgc.max_FeC_ratio_of_bacteria
    Bact = Bact(zₘₐₓ, z, Z, M) 
    Lₗᵢₘᵇᵃᶜᵗ = Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe, bgc)[2]
    return μₘₐₓ⁰*fₚ(T)*Lₗᵢₘᵇᵃᶜᵗ*θₘₐₓᶠᵉᵇᵃᶜᵗ*Fe*Bact/(K_Feᴮ¹ + Fe) 
end

@inline function (bgc::PISCES)(::Val{:Fe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust) #eq60
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    δᴾ = bgc.exudation_of_DOC.P
    δᴰ = bgc.exudation_of_DOC.D
    θᶠᵉᶻ = bgc.FeZ_redfield_ratio
    μₘₐₓ⁰ = bgc.growth_rate_at_zero
    θₘₐₓᶠᵉᴾ = bgc.max_iron_quota.P
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton.P
    K_Feᴾᶠᵉᵐⁱⁿ = bgc.min_half_saturation_const_for_iron_uptake.P
    Pₘₐₓ = bgc.threshold_concentration_for_size_dependency.P
    θₘₐₓᶠᵉᴰ = bgc.max_iron_quota.D
    Sᵣₐₜᴰ = bgc.size_ratio_of_phytoplankton.D
    K_Feᴰᶠᵉᵐⁱⁿ = bgc.min_half_saturation_const_for_iron_uptake.D
    Dₘₐₓ = bgc.threshold_concentration_for_size_dependency.D

    bFe = Fe

    L_Feᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[6]
    L_Feᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[6]
    
    sh = get_sh(z, zₘₓₗ)
    
    λₚₒ¹ = λ¹(T, O₂, bgc)

    μᴾᶠᵉ = μᴵᶠᵉ(P, Pᶠᵉ, θₘₐₓᶠᵉᴾ, Sᵣₐₜᴾ, K_Feᴾᶠᵉᵐⁱⁿ, Pₘₐₓ, L_Feᴾ, bFe, bgc)
    μᴰᶠᵉ = μᴵᶠᵉ(D, Dᶠᵉ, θₘₐₓᶠᵉᴰ, Sᵣₐₜᴰ, K_Feᴰᶠᵉᵐⁱⁿ, Dₘₐₓ, L_Feᴰ, bFe, bgc)
   
    #Iron quotas
    θᶠᵉᴾ = θ(Pᶠᵉ, P)
    θᶠᵉᴰ = θ(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = θ(SFe, POC)
    θᶠᵉᴳᴼᶜ = θ(BFe, GOC)
    #Grazing
    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazingᶻ(P, D, POC, T, bgc)
    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ = grazingᴹ(P, D, Z, POC, T, bgc)
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    
    ∑θᶠᵉⁱgᵢᶻ = θᶠᵉᴾ*grazingᶻ[2] + θᶠᵉᴰ*grazingᶻ[3] + θᶠᵉᴾᴼᶜ*grazingᶻ[4] #over P, D, POC
    ∑θᶠᵉⁱgᵢᴹ = θᶠᵉᴾ*grazingᴹ[2] + θᶠᵉᴰ*grazingᴹ[3] + θᶠᵉᴾᴼᶜ*grazingᴹ[4] + θᶠᵉᶻ*grazingᴹ[5] #graze on P, D, POC, Z 

    Bactfe = Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc)

    #Gross growth efficiency
    eₙᶻ = eₙᴶ(gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    eₙᴹ = eₙᴶ(gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ, Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC, bgc)
    
    return max(0, (1-σᶻ)*(∑θᶠᵉⁱgᵢᶻ/∑gᶻ - eₙᶻ*θᶠᵉᶻ))*∑gᶻ*Z 
    + max(0, (1-σᴹ)*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FF + θᶠᵉᴳᴼᶜ*g_GOC_FF )/(∑gᴹ+∑g_FFᴹ) - eₙᴹ*θᶠᵉᶻ)*(∑gᴹ+∑g_FFᴹ)*M 
    + γᴹ*θᶠᵉᶻ*Rᵤₚᴹ(M, T) + λₚₒ¹*SFe - (1 - δᴾ)*μᴾᶠᵉ*P - (1 - δᴰ)*μᴰᶠᵉ*D 
    - Scav(POC, GOC, CaCO₃, BSi, DOC, T, Fe) - Cgfe1(sh, Fe, POC, DOC, T, bgc) - Cgfe2(sh, Fe, T, DOC, GOC, bgc) 
    - Aggfe(Fe, DOC, T, bgc) - Bactfe
end

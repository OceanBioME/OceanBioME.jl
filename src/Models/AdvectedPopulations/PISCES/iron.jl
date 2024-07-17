
# TO DO: 
    #Fill in Bact function

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

@inline function Cgfe1(sh, Fe, POC, DOC, T)
    a₁ = bgc.aggregation_rate_of_DOC_to_POC_1
    a₂ = bgc.aggregation_rate_of_DOC_to_POC_2
    a₄ = bgc.aggregation_rate_of_DOC_to_POC_4
    a₅ = bgc.aggregation_rate_of_DOC_to_POC_5
   
    FeL = Fe - Fe¹(DOC, T, Fe) #eq64
    Fe_coll = 0.5*FeL
    return ((a₁*DOC + a₂*POC)*sh+a₄*POC + a₅*DOC)*Fe_coll
end

@inline function Cgfe2(sh, Fe, T, DOC, GOC)
    a₃ = bgc.aggregation_rate_of_DOC_to_GOC_3
    FeL = Fe - Fe¹(Fe, DOC, T)
    Fe_coll = 0.5*FeL
    return a₃*GOC*sh*Fe_coll
end

@inline function Aggfe(Fe, DOC, T)
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6)
    return 1000*λ_Fe*max(0, Fe - Lₜ)*Fe¹(DOC, T, Fe)
end

@inline function Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ)
    K_Feᴮ¹ = #add as parameter = 0.01 or 2.5e-10
    θₘₐₓᶠᵉᵇᵃᶜᵗ = #add as parameter, 10e-6?
    Bact = Bact(zₘₐₓ, z, Z, M) 
    Lₗᵢₘᵇᵃᶜᵗ = Lᵇᵃᶜᵗ(DOC, PO₄, NO₃, NH₄, bFe)[2]
    return μₘₐₓ⁰*fₚ(T)*Lₗᵢₘᵇᵃᶜᵗ*θₘₐₓᶠᵉᵇᵃᶜᵗ*Fe*Bact/(K_Feᴮ¹ + Fe) 
end

@inline function (pisces::PISCES)(::Val{:Fe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅) #eq60
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    δᴾ = bgc.exudation_of_DOC.P
    δᴰ = bgc.exudation_of_DOC.D
    θᶠᵉᶻ = bgc.FeZ_redfield_ratio
    μₘₐₓ⁰ = bgc.growth_rate_at_zero
    θₘₐₓᶠᵉᴾ = 
    Sᵣₐₜᴾ = 
    K_Feᴾᶠᵉᵐⁱⁿ = 
    Pₘₐₓ = 
    L_Feᴾ = 
    θₘₐₓᶠᵉᴰ = 
    Sᵣₐₜᴰ = 
    K_Feᴰᶠᵉᵐⁱⁿ = 
    Iₘₐₓ = 
    L_Feᴰ = 
    
    sh = get_sh(z, zₘₓₗ)
    
    λₚₒ¹ = λ¹(T, O₂)

    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC)

    μᴾᶠᵉ = μᴵᶠᵉ(P, Pᶠᵉ, θₘₐₓᶠᵉᴾ, Sᵣₐₜᴾ, K_Feᴾᶠᵉᵐⁱⁿ, Pₘₐₓ, L_Feᴾ, bFe)
    μᴰᶠᵉ = μᴵᶠᵉ(D, Dᶠᵉ, θₘₐₓᶠᵉᴰ, Sᵣₐₜᴰ, K_Feᴰᶠᵉᵐⁱⁿ, Iₘₐₓ, L_Feᴰ, bFe)
   
    #Iron quotas
    θᶠᵉᴾ = θ(Pᶠᵉ, P)
    θᶠᵉᴰ = θ(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = θ(SFe, POC)
    θᶠᵉᴳᴼᶜ = θ(BFe, GOC)
    #Grazing
    grazingᶻ = grazingᶻ(P, D, POC, T)
    grazingᴹ = grazingᴹ(P, D, Z, POC, T)
    ∑gᶻ = grazingᶻ[1]
    ∑gᴹ = grazingᴹ[1]
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)
    

    ∑θᶠᵉⁱgᵢᶻ = θᶠᵉᴾ*grazingᶻ[2] + θᶠᵉᴰ*grazingᶻ[3] + θᶠᵉᴾᴼᶜ*grazingᶻ[4] #over P, D, POC
    ∑θᶠᵉⁱgᵢᴹ = θᶠᵉᴾ*grazingᴹ[2] + θᶠᵉᴰ*grazingᴹ[3] + θᶠᵉᴾᴼᶜ*grazingᴹ[4] + θᶠᵉᶻ*grazingᴹ[5] #graze on P, D, POC, Z 

    Bactfe = Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ)

    return max(0, (1-σᶻ)*(∑θᶠᵉⁱgᵢᶻ/∑gᶻ - eₙᶻ*θᶠᵉᶻ)*∑gᶻ*Z + max(0, (1-σᴹ)*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FF + θᶠᵉᴳᴼᶜ*g_GOC_FF )/(∑gᴹ+∑g_FFᴹ) - eₙᴹ*θᶠᵉᶻ)*(∑gᴹ+∑g_FFᴹ)*M + γᴹ*θᶠᵉᶻ*Rᵤₚᴹ(M, T) + λₚₒ¹*SFe - (1 - δᴾ)*μᴾᶠᵉ*P - (1 - δᴰ)*μᴰᶠᵉ*D 
    - Scav(POC, GOC, CaCO₃, BSi, DOC, T, Fe) - Cgfe1(sh, Fe, POC, DOC, T) - Cgfe2(sh, Fe, T, DOC, GOC) - Aggfe(Fe, DOC, T) - Bactfe
end

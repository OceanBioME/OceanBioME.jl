
# TO DO: 
    #Should Lₜ be defined as in (67) or from the parameter list? Best way to define this for Aggfe, should we define internally in Fe¹()?
    #Where are θₘₐₓᶠᵉᵇᵃᶜᵗ?
    #Where is K_Feᴮ¹ defined? where is θₘₐₓᶠᵉᵇᵃᶜᵗ defined?

# Using simple chemistry model. 
# This document contains functions for the following:
    # Fe¹ (eq65), dissolved free inorganic iron
    # Cgfe1, Cgfe2, Aggfe, Bactfe (eqs 61, 62, 63)
    # Forcing for Fe (eq60)

@inline function Fe¹(DOC, T, Fe)
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6) # bgc.total_concentration_of_iron_ligands
    K_eqᶠᵉ = 10^(16.27 - 1565.7/max(T, 5))
    Δ = 1 +  K_eqᶠᵉ(T)*Lₜ -  K_eqᶠᵉ(T)*Fₑ

    return (-Δ + sqrt(Δ^2 + 4*K_eqᶠᵉ(T)*Fe))/2*K_eqᶠᵉ(T) #eq65
end

@inline function Cgfe1(DOC, POC, Fe, T)
    a₁ = bgc.aggregation_rate_of_DOC_to_POC_1
    a₂ = bgc.aggregation_rate_of_DOC_to_POC_2
    a₄ = bgc.aggregation_rate_of_DOC_to_POC_4
    a₅ = bgc.aggregation_rate_of_DOC_to_POC_5
   

    FeL = Fe - Fe¹(DOC, T, Fe)
    Fe_coll = 0.5*FeL
    return ((a₁*DOC + a₂*POC)*sh+a₄*POC + a₅*DOC)*Fe_coll
end

@inline function  Cgfe2(GOC, Fe, T, DOC)
    a₃ = bgc.aggregation_rate_of_DOC_to_GOC_3
    sh =

    FeL = Fe - Fe¹(DOC, T, Fe)
    Fe_coll = 0.5*FeL
    return a₃*GOC*sh*Fe_coll
end
@inline function Aggfe(Fe,DOC)
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    Lₜ = 
    return 1000*λ_Fe*max(0, Fe - Lₜ)*Fe¹(DOC, T, Fe)
end

@inline function (pisces::PISCES)(::Val{:Fe}, x, y, z, t, P, PAR) #(60)
    σᶻ = bgc.non_assimilated_fraction[1]
    γᴹ = bgc.excretion_as_DOM[2]
    σᴹ = bgc.non_assimilated_fraction[2]
    δᴾ = bgc.exudation_of_DOC[1]
    δᴰ = bgc.exudation_of_DOC[2]
   
    sh = # Find value
    θₘₐₓᶠᵉᵇᵃᶜᵗ = # check where is this defined?
    λₚₒ¹ = λ¹(T, O₂)

    Bact = 
    K_Feᴮ¹ = 

    Bactfe = μₚ()*Lₗᵢₘᵇᵃᶜᵗ()*θₘₐₓᶠᵉᵇᵃᶜᵗ*Fe*Bact/(K_Feᴮ¹ + Fe) # where is K_Feᴮ¹ defined? where is θₘₐₓᶠᵉᵇᵃᶜᵗ defined?

    return max(0, (1-σᶻ)*(∑θᶠᵉⁱ*gᶻ())/∑gᶻ() - eₙᶻ()θ(Zᶠᵉ, Z))*∑gᶻ()*Z + max(0, (1-σᴹ)*(∑θᶠᵉⁱ*gᴹ() + ∑θᶠᵉⁱ*g_FFᴹ())/(∑gᴹ()+g_FFᴹ()+g_FFᴹ()) - eₙᴹ()*θ(Zᶠᵉ, Z))*(∑gᴹ()+g_FFᴹ()+g_FFᴹ())*M + γᴹ*θ(Zᶠᵉ, Z)*Rᵤₚᴹ() + λ_poc1*SFe - (1 - δᴾ)*μᴾᶠᵉ()*P - (1 - δᴰ)*μᴰᶠᵉ()*D - Scav() - Cfge1() - Cgfe2() - Aggfe() - Bactfe()
end

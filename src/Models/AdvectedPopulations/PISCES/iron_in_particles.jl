
#TO DO:
    #How to define quotas, POCᶠᵉ is not defined?
    #Add partial derivative
    #D_dust?
    #Where is κ_Bactˢᶠᵉ, κ_Bactᴮᶠᵉ defined?
    #Change ω to w.
    #Where to split longer functions to make more readable?

#This document contains functions for the following:
    #θᶠᵉ for calculating iron quotas
    #Scav (eq50)
    #Forcing equations for SFe and BFe. (eqs 48 and 49)

@inline θᶠᵉ(J, Jᶠᵉ) = Jᶠᵉ/J #Is this the correct interpretation of the quotas? How does this work for POC?

@inline function λ_Fe¹(POC, GOC, CaCO₃, BSi) 
    λ_Feᵐⁱⁿ = bgc.min_scavenging_rate_of_iron
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    λ_Feᵈᵘˢᵗ = bgc.scavenging_rate_of_iron_by_dust
    ω_dust = bgc.sinking_speed_of_dust

    Dust = D_dust/ω_dust
    
    return λ_Feᵐⁱⁿ + λ_Fe*(POC, GOC, CaCO₃, BSi) + λ_Feᵈᵘˢᵗ*Dust #eq50
end

@inline Scav(POC, GOC, CaCO₃, BSi, DOC, T, Fe) = λ_Fe¹(POC, GOC, CaCO₃, BSi)*Fe¹(DOC, T, Fe)

@inline function (pisces::PISCES)(::Val{:SFe}, x, y, z, t, P, PAR) 
    
    σᶻ = bgc.non_assimilated_fraction[1]
    rᶻ = bgc.zooplankton_linear_mortality[1]
    mᶻ = bgc.zooplankton_quadratic_mortality[1]
    λ_GOC¹ = #where is this defined?
    mᴾ = bgc.zooplankton_quadratic_mortality[2]
    ωᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    mᴰ = bgc.phytoplankton_mortality_rate[2] 
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    λₚₒ¹ = λ¹(T, O₂)
    κ_Bactˢᶠᵉ = #where defined?
    ωₚₒ = bgc.sinking_speed_of_POC

    Fe¹ = Fe¹(DOC, T, Fe)

    grazingᶻ = grazingᶻ()
    ∑θᶠᵉⁱgᵢᶻ = θᶠᵉ(P, )*grazingᶻ[2] + θᶠᵉ(D, )*grazingᶻ[3] + θᶠᵉ(POC, )*grazingᶻ[4] #over P, D, POC

    g_POC_FFᴹ = g_FF*bₘ^T*w_POC*POC

    partial_CaCO₃ = 

    return σᶻ*∑θᶠᵉⁱgᵢᶻ*Z + θᶠᵉ(Z, )*(rᶻ*(b_Z^T)*K_mondo(Z, Kₘ)*Z + mᶻ*(b_Z^T)*(Z^2)) 
    + λ_GOC¹*BFe + θᶠᵉ(P, )*(1 - 0.5*R_CaCO₃())*(mᵖ*K_mondo(P, Kₘ)*P + sh*ωᴾ*P^2) + θᶠᵉ(D, )*0.5*mᴰ*K_mondo(D, Kₘ)*D + λ_Fe*POC*Fe¹ 
    + Cgfe1() - λₚₒ¹*SFe - θᶠᵉ(POC, )*Φ() - θᶠᵉ(POC, )*(grazingᴹ[4] + g_POC_FFᴹ)*M + κ_Bactˢᶠᵉ*Bactfe() - θᶠᵉ(POC, )*grazingᶻ()[4] - ωₚₒ*partial_CaCO₃ #add partial derivative #eq48
end 

@inline function (pisces::PISCES)(::Val{:BFe}, x, y, z, t, P, PAR) 
    
σᴹ = bgc.non_assimilated_fraction[2]
rᴹ = bgc.zooplankton_linear_mortality
mᴾ = bgc.phytoplankton_mortality_rate[1]
Kₘ = bgc.half_saturation_const_for_mortality
ωᴾ = bgc.min_quadratic_mortality_of_phytoplankton
mᴰ = bgc.phytoplankton_mortality_rate[2]
κ_Bactᴮᶠᵉ = #where defined?
λ_Fe = bgc.slope_of_scavenging_rate_of_iron
g_FF = bgc.flux_feeding_rate
ωₚₒ = bgc.sinking_speed_of_POC
bₘ = bgc.temperature_sensitivity_term[2]
ωₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms

Lₗᵢₘᴰ = 
ωᴰ = ωᴾ + ωₘₐₓᴰ*(1 - Lₗᵢₘᴰ)

Fe¹ = Fe¹(DOC, T, Fe)

grazingᴹ = grazingᴹ()
∑θᶠᵉⁱgᵢᴹ = θᶠᵉ(P, )*grazingᴹ[2] + θᶠᵉ(D, )*grazingᴹ[3] + θᶠᵉ(POC, )*grazingᴹ[4] + θᶠᵉ(Z )*grazingᴹ[5] #graze on P, D, POC, Z 

gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC 
g_GOC_FFᴹ = g_FF*bₘ^T*ω_GOC()*GOC 

    return σᴹ*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉ(POC,)*gₚₒ_FFᴹ + θᶠᵉ(GOC,)*g_GOC_FFᴹ)*M 
    + θᶠᵉ(M, )*(rᴹ*(bₘ^T)*K_mondo(M, Kₘ)*M + Pᵤₚᴹ()) + θᶠᵉ(P, )*0.5*R_CaCO₃()*(mᴾ*K_mondo(P, Kₘ)*P + sh*ωᴾ*P^2) 
    + θᶠᵉ(D, )*(0.5*mᴰ*K_mondo(D, Kₘ)*D + sh*ωᴰ*D^2) 
    + κ_Bactᴮᶠᵉ*Bactfe() + λ_Fe*GOC*Fe¹ + θᶠᵉ(POC, )*Φ() + Cgfe2() - θᶠᵉ(GOC, )* g_GOC_FFᴹ - λₚₒ¹*BFe - ω_GOC* #Add partial derivative
end
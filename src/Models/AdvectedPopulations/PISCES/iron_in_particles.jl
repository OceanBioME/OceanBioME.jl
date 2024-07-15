
#TO DO: Section unfinished - finish filling in symbols and check the equations
    #Check definition of fₘ(T)
    #Add partial derivative
    #What is Dust()?
    #Where is κ_Bactˢᶠᵉ, κ_Bactᴮᶠᵉ defined?
    #How to define quotas, POCᶠᵉ is not defined?
    #Where to split longer functions to make more readable?

#This document contains functions for the following:
    #θᶠᵉ for calculating iron quotas
    #Scav (eq50)
    #Forcing equations for SFe and BFe. (eqs 48 and 49)

@inline θᶠᵉ(J, Jᶠᵉ) = J/Jᶠᵉ #Is this the correct interpretation of the quotas? How does this work for POC?

@inline λ_Fe¹(POC, GOC, CaCO₃, BSi) 
    λ_Feᵐⁱⁿ = bgc.min_scavenging_rate_of_iron
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    λ_Feᵈᵘˢᵗ = bgc.scavenging_rate_of_iron_by_dust
    
    return  λ_Feᵐⁱⁿ + λ_Fe*(POC, GOC, CaCO₃, BSi) + λ_Feᵈᵘˢᵗ*Dust() #eq50, what is Dust()?
end

@inline Scav(POC, GOC, CaCO₃, BSi, DOC, T, Fe) = λ_Fe¹(POC, GOC, CaCO₃, BSi)*Fe¹(DOC, T, Fe)

@inline function (pisces::PISCES)(::Val{:SFe}, x, y, z, t, P, PAR) 
    
    σᶻ = 
    rᶻ = 
    mᶻ = 
    λ_GOC¹ = #where is this defined?
    mᴾ = 
    ωᴾ = 
    mᴰ = 
    λ_Fe = 
    λₚₒ¹ = 
    κ_Bactˢᶠᵉ = #where defined?
    ωₚₒ = 

    grazingᶻ = grazingᶻ()
    ∑θᶠᵉgᶻ = θᶠᵉ(P, )*grazingᶻ[2] + θᶠᵉ(D, )*grazingᶻ[3] + θᶠᵉ(POC, )*grazingᶻ[4] #over P, D, POC

    g_POC_FFᴹ = g_FF*bₘ^T*w_POC*POC

    return σᶻ*∑θᶠᵉgᶻ*Z + θᶠᵉ(Z, )*rᶻ*b_Z^T*Z^2 + λ_GOC¹*BFe + θᶠᵉ(P, )*(1 - 0.5*R_CaCO₃())*(mᵖ*K_mondo(P, Kₘ)*P + sh*ωᴾ*P^2) + θᶠᵉ(D, )*0.5*mᴰ*K_mondo(D, Kₘ)*D + λ_Fe*POC*Fe¹ + Cgfe1() - λ_POC¹*SFe - θᶠᵉ(POC, )*Φ() - θᶠᵉ(POC, )*(grazingᴹ[4] + g_POC_FFᴹ)*M + κ_Bactˢᶠᵉ*Bactfe() - θᶠᵉ(POC, )*grazingᶻ()[4] - ωₚₒ* #add partial derivative #eq48
end

@inline function (pisces::PISCES)(::Val{:BFe}, x, y, z, t, P, PAR) 
    
σᴹ = 
rᴹ = 
mᴾ = 
Kₘ = 
ωᴾ = 
mᴰ = 
ωᴰ = 
κ_Bactᴮᶠᵉ = #where defined?
λ_Fe = 
g_FF = 
ωₚₒ = 
bₘ = 

grazingᴹ = grazingᴹ()
∑θᶠᵉgᴹ = 

gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC 
g_GOC_FFᴹ = g_FF*bₘ^T*ω_GOC()*GOC 

    return σᴹ*(∑θᶠᵉgᴹ + θᶠᵉ(POC,)*gₚₒ_FFᴹ + θᶠᵉ(GOC,)*g_GOC_FFᴹ)*M + θᶠᵉ(M, )*(rᴹ*(bₘ^T)*K_mondo(M, Kₘ)*M + Pᵤₚᴹ()) + θᶠᵉ(P, )*0.5*R_CaCO₃()*(mᴾ*K_mondo(P, Kₘ)*P + sh*ωᴾ*P^2) + θᶠᵉ(D, )*(0.5*mᴰ*K_mondo(D, Kₘ)*D + sh*ωᴰ*D^2) + κ_Bactᴮᶠᵉ*Bactfe() + λ_Fe*GP*Fe¹ + θᶠᵉ(POC, )*Φ + Cgfe2() - θᶠᵉ(GOC, )* g_GOC_FFᴹ - λ_POC¹*BFe - ω_GOC* #Add partial derivative
end
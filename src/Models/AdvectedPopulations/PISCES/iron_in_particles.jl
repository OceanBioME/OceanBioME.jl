#This document contains functions for the following:
    #Scav (eq50)
    #Forcing equations for SFe and BFe. (eqs 48 and 49)


@inline function λ_Fe¹(POC, GOC, CaCO₃, BSi, D_dust) 
    λ_Feᵐⁱⁿ = bgc.min_scavenging_rate_of_iron
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    λ_Feᵈᵘˢᵗ = bgc.scavenging_rate_of_iron_by_dust
    w_dust = bgc.sinking_speed_of_dust

    Dust = D_dust/w_dust #eq84, check how to define D_dust?
    
    return λ_Feᵐⁱⁿ + λ_Fe*(POC + GOC + CaCO₃ + BSi) + λ_Feᵈᵘˢᵗ*Dust #eq50
end

@inline Scav(POC, GOC, CaCO₃, BSi, DOC, T, Fe) = λ_Fe¹(POC, GOC, CaCO₃, BSi)*Fe¹(DOC, T, Fe)

@inline function (pisces::PISCES)(::Val{:SFe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅) 
    #Parameters
    σᶻ = bgc.non_assimilated_fraction.Z
    rᶻ = bgc.zooplankton_linear_mortality.Z
    Kₘ = bgc.half_saturation_const_for_mortality
    mᶻ = bgc.zooplankton_quadratic_mortality.Z
    mᴾ = bgc.zooplankton_quadratic_mortality.P
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    mᴰ = bgc.phytoplankton_mortality_rate.D
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    κ_Bactˢᶠᵉ = bgc.coefficient_of_bacterial_uptake_of_iron_in_POC
    wₚₒ = bgc.sinking_speed_of_POC
    g_FF = bgc.flux_feeding_rate
    bₘ = bgc.temperature_sensitivity_term.M

    sh = get_sh(z, zₘₓₗ)

    Fe¹ = Fe¹(DOC, T, Fe)
    λₚₒ¹ = λ¹(T, O₂)
    #Iron quotas
    θᶠᵉᴾ = θ(Pᶠᵉ, P)
    θᶠᵉᴰ = θ(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = θ(SFe, POC)
    #Grazing
    grazingᶻ = grazingᶻ(P, D, POC, T)
    ∑θᶠᵉⁱgᵢᶻ = θᶠᵉᴾ*grazingᶻ[2] + θᶠᵉᴰ*grazingᶻ[3] + θᶠᵉᴾᴼᶜ*grazingᶻ[4] #over P, D, POC
    gₚₒ_FFᴹ = g_FF*(bₘ^T)*wₚₒ*POC
    #Bacteria iron
    Bactfe = Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ)

    return σᶻ*∑θᶠᵉⁱgᵢᶻ*Z + θᶠᵉᶻ*(rᶻ*(b_Z^T)*K_mondo(Z, Kₘ)*Z + mᶻ*(b_Z^T)*(Z^2)) + λₚₒ¹*BFe + θᶠᵉᴾ*(1 - 0.5*R_CaCO₃(P, T, PAR, zₘₓₗ))*(mᵖ*K_mondo(P, Kₘ)*P + sh*wᴾ*P^2) + θᶠᵉᴰ*0.5*mᴰ*K_mondo(D, Kₘ)*D + λ_Fe*POC*Fe¹ + Cgfe1(h, Fe, POC, DOC, T) - λₚₒ¹*SFe - θᶠᵉᴾᴼᶜ*Φ(POC, GOC, sh) - θᶠᵉᴾᴼᶜ*(grazingᴹ[4] + gₚₒ_FFᴹ)*M + κ_Bactˢᶠᵉ*Bactfe - θᶠᵉᴾᴼᶜ*grazingᶻ[4] #Partial derivative omitted #eq48
end 

@inline function (pisces::PISCES)(::Val{:BFe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅) 
    #Parameters
    σᴹ = bgc.non_assimilated_fraction.M
    rᴹ = bgc.zooplankton_linear_mortality
    mᴾ = bgc.phytoplankton_mortality_rate.P
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    mᴰ = bgc.phytoplankton_mortality_rate.D
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    g_FF = bgc.flux_feeding_rate
    wₚₒ = bgc.sinking_speed_of_POC
    bₘ = bgc.temperature_sensitivity_term.M
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms
    κ_Bactᴮᶠᵉ = bgc.coefficient_of_bacterial_uptake_of_iron_in_GOC
    w_GOCᵐⁱⁿ = bgc.min_sinking_speed_of_GOC

    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅)[1]

    wᴰ = wᴾ + wₘₐₓᴰ*(1 - Lₗᵢₘᴰ)

    Fe¹ = Fe¹(DOC, T, Fe)
    #Iron quotas
    θᶠᵉᴾ = θ(Pᶠᵉ, P)
    θᶠᵉᴰ = θ(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = θ(SFe, POC)
    θᶠᵉᴳᴼᶜ = θ(BFe, GOC)
    #Grazing
    grazingᴹ = grazingᴹ(P, D, Z, POC, T)
    ∑θᶠᵉⁱgᵢᴹ = θᶠᵉᴾ*grazingᴹ[2] + θᶠᵉᴰ*grazingᴹ[3] + θᶠᵉᴾᴼᶜ*grazingᴹ[4] + θᶠᵉᶻ*grazingᴹ[5] #graze on P, D, POC, Z 
    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC 
    zₘₐₓ = max(zₑᵤ, zₘₓₗ)   #41a
    w_GOC = w_GOCᵐⁱⁿ + (200 - w_GOCᵐⁱⁿ)*(max(0, z-zₘₐₓ))/(5000) #41b
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC 

    return σᴹ*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FFᴹ + θᶠᵉᴳᴼᶜ*g_GOC_FFᴹ)*M + θᶠᵉᴹ*(rᴹ*(bₘ^T)*K_mondo(M, Kₘ)*M + Pᵤₚᴹ()) + θᶠᵉᴾ*0.5*R_CaCO₃()*(mᴾ*K_mondo(P, Kₘ)*P + sh*wᴾ*P^2) + θᶠᵉᴰ*(0.5*mᴰ*K_mondo(D, Kₘ)*D + sh*wᴰ*D^2) + κ_Bactᴮᶠᵉ*Bactfe() + λ_Fe*GOC*Fe¹ + θᶠᵉᴾᴼᶜ*Φ() + Cgfe2() - θᶠᵉᴳᴼᶜ* g_GOC_FFᴹ*M - λₚₒ¹*BFe #Partial derivative omitted
end
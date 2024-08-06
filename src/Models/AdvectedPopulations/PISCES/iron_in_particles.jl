#This document contains functions for the following:
    #Scav (eq50)
    #Forcing equations for SFe and BFe. (eqs 48 and 49)


@inline function λ_Fe¹(POC, GOC, CaCO₃, PSi, D_dust, bgc) 
    λ_Feᵐⁱⁿ = bgc.min_scavenging_rate_of_iron
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    λ_Feᵈᵘˢᵗ = bgc.scavenging_rate_of_iron_by_dust
    w_dust = bgc.sinking_speed_of_dust

    Dust = D_dust/(w_dust+ eps(0.0)) #eq84, check how to define D_dust?
    
    return λ_Feᵐⁱⁿ + λ_Fe*(POC + GOC + CaCO₃ + PSi) + λ_Feᵈᵘˢᵗ*Dust #eq50
end

@inline Scav(POC, GOC, CaCO₃, PSi, D_dust, DOC, T, Fe, bgc) = λ_Fe¹(POC, GOC, CaCO₃, PSi, D_dust, bgc)*get_Fe¹(Fe, DOC, T)

@inline function (bgc::PISCES)(::Val{:SFe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³) 
    #Parameters
    σᶻ = bgc.non_assimilated_fraction.Z
    rᶻ = bgc.zooplankton_linear_mortality.Z
    Kₘ = bgc.half_saturation_const_for_mortality
    mᶻ = bgc.zooplankton_quadratic_mortality.Z
    mᴾ = bgc.phytoplankton_mortality_rate.P
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    mᴰ = bgc.phytoplankton_mortality_rate.D
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    κ_Bactˢᶠᵉ = bgc.coefficient_of_bacterial_uptake_of_iron_in_POC
    wₚₒ = bgc.sinking_speed_of_POC
    g_FF = bgc.flux_feeding_rate
    bₘ = bgc.temperature_sensitivity_term.M
    b_Z = bgc.temperature_sensitivity_term.Z
    μₘₐₓ⁰ = bgc.growth_rate_at_zero
    θᶠᵉᶻ = bgc.FeC_ratio_of_zooplankton

    sh = get_sh(z, zₘₓₗ)

    bFe = Fe
    zₘₐₓ = max(zₑᵤ, zₘₓₗ)
    Fe¹ = get_Fe¹(Fe, DOC, T)   #same name
    λₚₒ¹ = λ¹(T, O₂, bgc)
    #Iron quotas
    θᶠᵉᴾ = θ(Pᶠᵉ, P)
    θᶠᵉᴰ = θ(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = θ(SFe, POC)
    #Grazing
    grazingᶻ = get_grazingᶻ(P, D, POC, T, bgc)
    grazingᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)
    ∑θᶠᵉⁱgᵢᶻ = θᶠᵉᴾ*grazingᶻ[2] + θᶠᵉᴰ*grazingᶻ[3] + θᶠᵉᴾᴼᶜ*grazingᶻ[4] #over P, D, POC
    gₚₒ_FFᴹ = g_FF*(bₘ^T)*wₚₒ*POC
    #Bacteria iron
    Bactfe = get_Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc)

    println("Sum of positive terms in SFe is ", σᶻ*∑θᶠᵉⁱgᵢᶻ*Z + θᶠᵉᶻ*(rᶻ*(b_Z^T)*K_mondo(Z, Kₘ)*Z + mᶻ*(b_Z^T)*(Z^2)) + λₚₒ¹*BFe + θᶠᵉᴾ*(1 - 0.5*get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc))*(mᴾ*K_mondo(P, Kₘ)*P + sh*wᴾ*P^2) + θᶠᵉᴰ*0.5*mᴰ*K_mondo(D, Kₘ)*D + λ_Fe*POC*Fe¹ + Cgfe1(sh, Fe, POC, DOC, T, bgc) + κ_Bactˢᶠᵉ*Bactfe)
    println("term a = $(σᶻ*∑θᶠᵉⁱgᵢᶻ*Z)), term b = $( θᶠᵉᶻ*(rᶻ*(b_Z^T)*K_mondo(Z, Kₘ)*Z + mᶻ*(b_Z^T)*(Z^2))), term c = $( λₚₒ¹*BFe), term d = $(λₚₒ¹*BFe + θᶠᵉᴾ*(1 - 0.5*get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc))*(mᴾ*K_mondo(P, Kₘ)*P + sh*wᴾ*P^2)), term e = $(θᶠᵉᴰ*0.5*mᴰ*K_mondo(D, Kₘ)*D), term f = $(λ_Fe*POC*Fe¹), term g = $(Cgfe1(sh, Fe, POC, DOC, T, bgc)), term h = $(κ_Bactˢᶠᵉ*Bactfe)")
    println("Sum of negative terms is ", λₚₒ¹*SFe + θᶠᵉᴾᴼᶜ*get_Φ(POC, GOC, sh, bgc) + θᶠᵉᴾᴼᶜ*(grazingᴹ[4] + gₚₒ_FFᴹ)*M+ θᶠᵉᴾᴼᶜ*grazingᶻ[4])
    println("term 1 = $(λₚₒ¹*SFe), term 2 = $(θᶠᵉᴾᴼᶜ*get_Φ(POC, GOC, sh, bgc)), term 3 = $(θᶠᵉᴾᴼᶜ*(grazingᴹ[4] + gₚₒ_FFᴹ)*M), term 4 = $(θᶠᵉᴾᴼᶜ*grazingᶻ[4]))")
    println("--")
    println("Total change is ", σᶻ*∑θᶠᵉⁱgᵢᶻ*Z + θᶠᵉᶻ*(rᶻ*(b_Z^T)*K_mondo(Z, Kₘ)*Z + mᶻ*(b_Z^T)*(Z^2)) + λₚₒ¹*BFe + θᶠᵉᴾ*(1 - 0.5*get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc))*(mᴾ*K_mondo(P, Kₘ)*P + sh*wᴾ*P^2) + θᶠᵉᴰ*0.5*mᴰ*K_mondo(D, Kₘ)*D + λ_Fe*POC*Fe¹ + Cgfe1(sh, Fe, POC, DOC, T, bgc) - λₚₒ¹*SFe - θᶠᵉᴾᴼᶜ*get_Φ(POC, GOC, sh, bgc) - θᶠᵉᴾᴼᶜ*(grazingᴹ[4] + gₚₒ_FFᴹ)*M + κ_Bactˢᶠᵉ*Bactfe - θᶠᵉᴾᴼᶜ*grazingᶻ[4])
    println("-------------------------------------")

    return σᶻ*∑θᶠᵉⁱgᵢᶻ*Z + θᶠᵉᶻ*(rᶻ*(b_Z^T)*K_mondo(Z, Kₘ)*Z + mᶻ*(b_Z^T)*(Z^2)) + λₚₒ¹*BFe + θᶠᵉᴾ*(1 - 0.5*get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc))*(mᴾ*K_mondo(P, Kₘ)*P + sh*wᴾ*P^2) + θᶠᵉᴰ*0.5*mᴰ*K_mondo(D, Kₘ)*D + λ_Fe*POC*Fe¹ + Cgfe1(sh, Fe, POC, DOC, T, bgc) - λₚₒ¹*SFe - θᶠᵉᴾᴼᶜ*get_Φ(POC, GOC, sh, bgc) - θᶠᵉᴾᴼᶜ*(grazingᴹ[4] + gₚₒ_FFᴹ)*M + κ_Bactˢᶠᵉ*Bactfe - θᶠᵉᴾᴼᶜ*grazingᶻ[4]*Z #Partial derivative omitted #eq48
end 

@inline function (bgc::PISCES)(::Val{:BFe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³)
    #Parameters
    σᴹ = bgc.non_assimilated_fraction.M
    rᴹ = bgc.zooplankton_linear_mortality.M
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
    θᶠᵉᶻ = bgc.FeC_ratio_of_zooplankton
    μₘₐₓ⁰ = bgc.growth_rate_at_zero

    bFe = Fe

    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]

    wᴰ = wᴾ + wₘₐₓᴰ*(1 - Lₗᵢₘᴰ)

    Fe¹ = get_Fe¹(Fe, DOC, T)
    #Iron quotas
    θᶠᵉᴾ = θ(Pᶠᵉ, P)
    θᶠᵉᴰ = θ(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = θ(SFe, POC)
    θᶠᵉᴳᴼᶜ = θ(BFe, GOC)

    λₚₒ¹ =  λ¹(T, O₂, bgc)

    sh = get_sh(z, zₘₓₗ)
    #Grazing
    grazingᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)
    ∑θᶠᵉⁱgᵢᴹ = θᶠᵉᴾ*grazingᴹ[2] + θᶠᵉᴰ*grazingᴹ[3] + θᶠᵉᴾᴼᶜ*grazingᴹ[4] + θᶠᵉᶻ*grazingᴹ[5] #graze on P, D, POC, Z 
    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC 
    zₘₐₓ = max(zₑᵤ, zₘₓₗ)   #41a
    w_GOC = w_GOCᵐⁱⁿ + (200 - w_GOCᵐⁱⁿ)*(max(0, z-zₘₐₓ))/(5000) #41b
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC 

    return σᴹ*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FFᴹ + θᶠᵉᴳᴼᶜ*g_GOC_FFᴹ)*M + θᶠᵉᶻ*(rᴹ*(bₘ^T)*K_mondo(M, Kₘ)*M + Pᵤₚ(M, T, bgc)) + θᶠᵉᴾ*0.5*get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc)*(mᴾ*K_mondo(P, Kₘ)*P + sh*wᴾ*P^2) + θᶠᵉᴰ*(0.5*mᴰ*K_mondo(D, Kₘ)*D + sh*wᴰ*D^2) + κ_Bactᴮᶠᵉ*get_Bactfe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc) + λ_Fe*GOC*Fe¹ + θᶠᵉᴾᴼᶜ*get_Φ(POC, GOC, sh, bgc) + Cgfe2(sh, Fe, T, DOC, GOC, bgc) - θᶠᵉᴳᴼᶜ* g_GOC_FFᴹ*M - λₚₒ¹*BFe #Partial derivative omitted
end
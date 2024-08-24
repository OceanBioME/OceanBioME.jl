#This document contains functions for the following:
    #Scav (eq50)
    #Forcing equations for SFe and BFe. (eqs 48 and 49)
#We use the 2 compartment version of the model, with iron comparments for small and big particles.

#Free form of iron is the only form susceptible to scavenging. We formulate the scavenging rate.
#Iron is scavenged by lithogenic and biogenic particles. Iron scavenged by POC and GOC routed to SFe and BFe, al other scavenging lost from the system.
@inline function Fe_scavenging_rate(POC, GOC, CaCO₃, PSi, D_dust, bgc) 
    λ_Feᵐⁱⁿ = bgc.min_scavenging_rate_of_iron
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    λ_Feᵈᵘˢᵗ = bgc.scavenging_rate_of_iron_by_dust
    w_dust = bgc.sinking_speed_of_dust
    Dust = D_dust/(w_dust+ eps(0.0)) #eq84
    return λ_Feᵐⁱⁿ + λ_Fe*(POC + GOC + CaCO₃ + PSi) + λ_Feᵈᵘˢᵗ*Dust #eq50
end

#Scavenging of free form of dissolved iron.
@inline Fe_scavenging(POC, GOC, CaCO₃, PSi, D_dust, DOC, T, Fe, bgc) = Fe_scavenging_rate(POC, GOC, CaCO₃, PSi, D_dust, bgc)*free_organic_iron(Fe, DOC, T)

@inline function (bgc::PISCES)(::Val{:SFe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃) 
    #Parameters
    σᶻ = bgc.non_assimilated_fraction.Z
    rᶻ = bgc.zooplankton_linear_mortality.Z
    Kₘ = bgc.half_saturation_const_for_mortality
    mᴾ, mᴰ = bgc.phytoplankton_mortality_rate
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    mᶻ = bgc.zooplankton_quadratic_mortality.Z
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    κ_Bactˢᶠᵉ = bgc.coefficient_of_bacterial_uptake_of_iron_in_POC
    wₚₒ = bgc.sinking_speed_of_POC
    g_FF = bgc.flux_feeding_rate
    b_Z, bₘ = bgc.temperature_sensitivity_term
    μₘₐₓ⁰ = bgc.growth_rate_at_zero
    θᶠᵉᶻ = bgc.FeC_ratio_of_zooplankton

    #Also required
    sh = shear_rate(z, zₘₓₗ)
    Fe¹ = free_organic_iron(Fe, DOC, T) 
    λₚₒ¹ = particles_carbon_degradation_rate(T, O₂, bgc)
    bFe = Fe
    
    #Iron quotas
    θᶠᵉᴾ = nutrient_quota(Pᶠᵉ, P)
    θᶠᵉᴰ = nutrient_quota(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = nutrient_quota(SFe, POC)
    
    #Grazing
    grazingᶻ = grazing_Z(P, D, POC, T, bgc)
    grazingᴹ = grazing_M(P, D, Z, POC, T, bgc)
    ∑θᶠᵉⁱgᵢᶻ = θᶠᵉᴾ*grazingᶻ[2] + θᶠᵉᴰ*grazingᶻ[3] + θᶠᵉᴾᴼᶜ*grazingᶻ[4] #over P, D, POC
    gₚₒ_FFᴹ = g_FF*(bₘ^T)*wₚₒ*POC
    
    #Bacteria iron
    zₘₐₓ = max(abs(zₑᵤ), abs(zₘₓₗ))
    Bactfe = bacterial_uptake_Fe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc)

    return (σᶻ*∑θᶠᵉⁱgᵢᶻ*Z 
            + θᶠᵉᶻ*(rᶻ*(b_Z^T)*(concentration_limitation(Z, Kₘ) + 3*oxygen_conditions(O₂, bgc))*Z + mᶻ*(b_Z^T)*(Z^2)) 
            + λₚₒ¹*BFe 
            + θᶠᵉᴾ*(1 - 0.5*rain_ratio(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc))*(mᴾ*concentration_limitation(P, Kₘ)*P + sh*wᴾ*P^2) 
            + θᶠᵉᴰ*0.5*mᴰ*concentration_limitation(D, Kₘ)*D + λ_Fe*POC*Fe¹ 
            + iron_colloid_aggregation_1(sh, Fe, POC, DOC, T, bgc) - λₚₒ¹*SFe - θᶠᵉᴾᴼᶜ*POC_aggregation(POC, GOC, sh, bgc) 
            - θᶠᵉᴾᴼᶜ*(grazingᴹ[4] + gₚₒ_FFᴹ)*M 
            + κ_Bactˢᶠᵉ*Bactfe - θᶠᵉᴾᴼᶜ*grazingᶻ[4]*Z) #eq48, partial derivative ommitted

    #Changes made from paper:
        #3*oxygen_conditions added to zooplankton linear mortality terms.
        #Z factor missing from final term.

end 

@inline function (bgc::PISCES)(::Val{:BFe}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    σᴹ = bgc.non_assimilated_fraction.M
    rᴹ = bgc.zooplankton_linear_mortality.M
    mᴾ, mᴰ = bgc.phytoplankton_mortality_rate
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron
    g_FF = bgc.flux_feeding_rate
    wₚₒ = bgc.sinking_speed_of_POC
    bₘ = bgc.temperature_sensitivity_term.M
    κ_Bactᴮᶠᵉ = bgc.coefficient_of_bacterial_uptake_of_iron_in_GOC
    θᶠᵉᶻ = bgc.FeC_ratio_of_zooplankton
    μₘₐₓ⁰ = bgc.growth_rate_at_zero

    #Other required terms
    sh = shear_rate(z, zₘₓₗ)
    Fe¹ = free_organic_iron(Fe, DOC, T) 
    λₚₒ¹ = particles_carbon_degradation_rate(T, O₂, bgc)
    wᴰ = D_quadratic_mortality(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)
    bFe = Fe
    
    #Iron quotas
    θᶠᵉᴾ = nutrient_quota(Pᶠᵉ, P)
    θᶠᵉᴰ = nutrient_quota(Dᶠᵉ, D)
    θᶠᵉᴾᴼᶜ = nutrient_quota(SFe, POC)
    θᶠᵉᴳᴼᶜ = nutrient_quota(BFe, GOC)

    #Grazing
    grazingᴹ = grazing_M(P, D, Z, POC, T, bgc)
    ∑θᶠᵉⁱgᵢᴹ = θᶠᵉᴾ*grazingᴹ[2] + θᶠᵉᴰ*grazingᴹ[3] + θᶠᵉᴾᴼᶜ*grazingᴹ[4] + θᶠᵉᶻ*grazingᴹ[5] #graze on P, D, POC, Z 
    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC 
    zₘₐₓ = max(abs(zₑᵤ), abs(zₘₓₗ))   #41a
    w_GOC = sinking_speed_of_GOC(z, zₑᵤ, zₘₓₗ, bgc)
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC 

    return (σᴹ*(∑θᶠᵉⁱgᵢᴹ + θᶠᵉᴾᴼᶜ*gₚₒ_FFᴹ + θᶠᵉᴳᴼᶜ*g_GOC_FFᴹ)*M 
            + θᶠᵉᶻ*(rᴹ*(bₘ^T)*(concentration_limitation(M, Kₘ) + 3*oxygen_conditions(O₂, bgc))*M + production_of_fecal_pellets(M, T, bgc)) 
            + θᶠᵉᴾ*0.5*rain_ratio(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc)*(mᴾ*concentration_limitation(P, Kₘ)*P + sh*wᴾ*P^2) 
            + θᶠᵉᴰ*(0.5*mᴰ*concentration_limitation(D, Kₘ)*D + sh*wᴰ*D^2) 
            + κ_Bactᴮᶠᵉ*bacterial_uptake_Fe(μₘₐₓ⁰, z, Z, M, Fe, DOC, PO₄, NO₃, NH₄, bFe, T, zₘₐₓ, bgc) 
            + λ_Fe*GOC*Fe¹ + θᶠᵉᴾᴼᶜ*POC_aggregation(POC, GOC, sh, bgc) + iron_colloid_aggregation_2(sh, Fe, T, DOC, GOC, bgc) 
            - θᶠᵉᴳᴼᶜ* g_GOC_FFᴹ*M - λₚₒ¹*BFe) #eq49, partial derivative omitted

    #Changes made from paper:
        #3*oxygen_conditions added to zooplankton linear mortality terms.
end
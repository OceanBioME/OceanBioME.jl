
#Particles of carbon are significant as sink and export carbon to the deep ocean.
#GOC are faster sinking particles with variable sinking speed.

# This documeent contains functions for:
    # Φ (eq39)
    # POC, GOC (eqs37, 40)

#Aggregation of POC due to turbulence and differential settling
@inline function POC_aggregation(POC, GOC, sh, bgc)
    a₆ = bgc.aggregation_rate_of_POC_to_GOC_6
    a₇ = bgc.aggregation_rate_of_POC_to_GOC_7
    a₈ = bgc.aggregation_rate_of_POC_to_GOC_8
    a₉ = bgc.aggregation_rate_of_POC_to_GOC_9

    return sh*a₆*POC^2 + sh*a₇*POC*GOC + a₈*POC*GOC + a₉*POC^2 #eq39
end

#Forcing for POC
@inline function (bgc::PISCES)(::Val{:POC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    σᶻ = bgc.non_assimilated_fraction.Z
    mᴾ, mᴰ = bgc.phytoplankton_mortality_rate
    mᶻ = bgc.zooplankton_quadratic_mortality.Z
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    wₚₒ = bgc.sinking_speed_of_POC
    rᶻ = bgc.zooplankton_linear_mortality.Z
    Kₘ = bgc.half_saturation_const_for_mortality
    b_Z, bₘ = bgc.temperature_sensitivity_term
    g_FF = bgc.flux_feeding_rate

    #Grazing
    grazing = grazing_Z(P, D, POC, T, bgc)
    ∑gᶻ = grazing[1]
    gₚₒᶻ = grazing[4]
    gₚₒᴹ = grazing_M(P, D, Z, POC, T, bgc)[4]
    gₚₒ_FFᴹ = g_FF*(bₘ^T)*wₚₒ*POC #29a
    
    #Aggregation
    sh = shear_rate(z, zₘₓₗ)
    Φ₁ᴰᴼᶜ = aggregation_process_for_DOC(DOC, POC, GOC, sh, bgc)[1]
    Φ₃ᴰᴼᶜ = aggregation_process_for_DOC(DOC, POC, GOC, sh, bgc)[3]
    Φ = POC_aggregation(POC, GOC, sh, bgc)

    R_CaCO₃ = rain_ratio(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc) 
    λₚₒ¹ = particles_carbon_degradation_rate(T, O₂, bgc)

    return (σᶻ*∑gᶻ*Z + 0.5*mᴰ*concentration_limitation(D, Kₘ)*D + rᶻ*(b_Z^T)*(concentration_limitation(Z, Kₘ)
            + 3*oxygen_conditions(O₂, bgc))*Z + mᶻ*(b_Z^T)*Z^2 
            + (1 - 0.5*R_CaCO₃)*(mᴾ*concentration_limitation(P, Kₘ)*P + sh*wᴾ*P^2) 
            + λₚₒ¹*GOC + Φ₁ᴰᴼᶜ + Φ₃ᴰᴼᶜ - (gₚₒᴹ + gₚₒ_FFᴹ)*M - gₚₒᶻ*Z - λₚₒ¹*POC - Φ)  #eq37, partial derivative ommitted as included elsewhere in OceanBioME
end

#Forcing for GOC
@inline function (bgc::PISCES)(::Val{:GOC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃)
    #Parameters
    σᴹ = bgc.non_assimilated_fraction.M
    mᴾ, mᴰ = bgc.phytoplankton_mortality_rate
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    rᴹ = bgc.zooplankton_linear_mortality.M
    Kₘ = bgc.half_saturation_const_for_mortality
    bₘ = bgc.temperature_sensitivity_term.M
    g_FF = bgc.flux_feeding_rate
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms
    wₚₒ = bgc.sinking_speed_of_POC

    #Grazing
    w_GOC = sinking_speed_of_GOC(z, zₑᵤ, zₘₓₗ, bgc)
    ∑gᴹ = grazing_M(P, D, Z, POC, T, bgc)[1] 
    ∑g_FFᴹ, gₚₒ_FFᴹ, g_GOC_FFᴹ = flux_feeding(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    
    #Aggregation
    sh = shear_rate(z, zₘₓₗ)
    Φ = POC_aggregation(POC, GOC, sh, bgc)
    Φ₂ᴰᴼᶜ = aggregation_process_for_DOC(DOC, POC, GOC, sh, bgc)[2]

    Pᵤₚᴹ = production_of_fecal_pellets(M, T, bgc)
    R_CaCO₃ = rain_ratio(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc)
    λₚₒ¹ = particles_carbon_degradation_rate(T, O₂, bgc)
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    wᴰ =  D_quadratic_mortality(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)

    return (σᴹ*(∑gᴹ + ∑g_FFᴹ)*M + rᴹ*(bₘ^T)*(concentration_limitation(M, Kₘ) + 3*oxygen_conditions(O₂, bgc))*M 
            + Pᵤₚᴹ + 0.5*R_CaCO₃*(mᴾ*concentration_limitation(P, Kₘ)*P + sh*wᴾ*P^2) + 0.5*mᴰ*concentration_limitation(D, Kₘ)*D
            + sh*D^2*wᴰ + Φ + Φ₂ᴰᴼᶜ - g_GOC_FFᴹ*M - λₚₒ¹*GOC) #eq40, partial derivative ommitted as included elsewhere in OceanBioME
end
# Shear rate still to be added
# Derivates also to still be added

@inline function Φ(POC, GOC, sh)
    a₆ = bgc.aggregation_rate_of_POC_to_GOC_6
    a₇ = bgc.aggregation_rate_of_POC_to_GOC_7
    a₈ = bgc.aggregation_rate_of_POC_to_GOC_8
    a₉ = bgc.aggregation_rate_of_POC_to_GOC_9

    return sh*a₆*POC^2 + sh*a₇*POC*GOC + a₈*POC*GOC + a₉*POC^2  #39
end

@inline function (pisces::PISCES)(::Val{:POC}, x, y, z, t, DOC, P, D, Pᶜʰˡ, Dᶜʰˡ, N, Fe, O₂, NO₃, PARᴾ, PARᴰ, Z, M, POC, GOC, T, L_day, zₘₓₗ, zₑᵤ)
    σᶻ = bgc.non_assimilated_fraction[1]
    mᴾ = bgc.phytoplankton_mortality_rate[1]
    mᶻ = bgc.zooplankton_quadratic_mortality[1]
    mᴰ = bgc.phytoplankton_mortality_rate[2]
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    wₚₒ = bgc.sinking_speed_of_POC
    rᶻ = bgc.zooplankton_linear_mortality[1]
    Kₘ = bgc.half_saturation_const_for_mortality
    b_z, bₘ = bgc.temperature_sensitivity_term
    g_FF = bgc.flux_feeding_rate

    grazing = grazingᶻ(P, D, POC, T)
    ∑gᶻ = grazing[1]
    gₚₒᶻ = grazing[4]

    sh =

    R_CaCO₃ = R_CaCO₃(zₘₓₗ, T, P, PAR) 
    λₚₒ¹ = λ¹(T, O₂)
    Φ₁ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh)[1]
    Φ₃ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh)[3]

    gₚₒᴹ = grazingᴹ(P, D, Z, POC, T)[4]
    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC #29a
    Φ = Φ(POC, GOC, sh)

    dPOCdz = 

    return σᶻ*∑gᶻ*Z + 0.5*mᴰ*K_mondo(D, Kₘ) + rᶻ*b_z^T*K_mondo(Z, Kₘ)*Z +
     mᶻ*b_z^T*Z^2 + (1 - 0.5*R_CaCO₃)*(mᴾ*K_mondo(P, Kₘ)*P + wᴾ*P^2) + 
     λₚₒ¹*GOC + Φ₁ᴰᴼᶜ + Φ₃ᴰᴼᶜ - (gₚₒᴹ + gₚₒ_FFᴹ)*M - gₚₒᶻ*Z - λₚₒ¹*POC - Φ - wₚₒ*dPOCdz  #37
end

@inline function (pisces::PISCES)(::Val{:GOC}, x, y, z, t, DOC, P, D, Pᶜʰˡ, Dᶜʰˡ, N, Fe, O₂, NO₃, PARᴾ, PARᴰ, Z, M, POC, GOC, T, L_day, zₘₓₗ, zₑᵤ)
    σᴹ = bgc.non_assimilated_fraction[2]
    mᴾ = bgc.phytoplankton_mortality_rate[1]
    mᴰ = bgc.phytoplankton_mortality_rate[2]
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    rᴹ = bgc.zooplankton_linear_mortality[2]
    Kₘ = bgc.half_saturation_const_for_mortality
    bₘ = bgc.temperature_sensitivity_term[2]
    g_FF = bgc.flux_feeding_rate
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms

    ∑gᴹ = grazingᴹ(P, D, Z, POC, T)[1] 
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)
    
    Pᵤₚᴹ = Pᵤₚ(M, T)
    R_CaCO₃ = R_CaCO₃(zₘₓₗ, T, P, PAR)

    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[1]
    wᴰ =  wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ)
    Φ = Φ(POC, GOC, sh)
    Φ₂ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh)[2]
    w_GOC = w(zₑᵤ, zₘₓₗ)
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC #29b
    λₚₒ¹ = λ¹(T, O₂)
    
    dGOCdz =

    return σᴹ*(∑gᴹ + ∑g_FFᴹ)*M + rᴹ*bₘ^T*K_mondo(M, Kₘ)*M + Pᵤₚᴹ + 
    0.5*R_CaCO₃*(mᴾ*K_mondo(P, Kₘ)*P + wᴾ*P^2) + 0.5*mᴰ*K_mondo(D, Kₘ)*D^3*wᴰ +
     Φ + Φ₂ᴰᴼᶜ - g_GOC_FFᴹ*M - λₚₒ¹*GOC - w_GOC*dGOCdz    #40
end
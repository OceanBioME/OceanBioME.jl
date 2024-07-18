@inline function Φ(POC, GOC, sh, bgc)
    a₆ = bgc.aggregation_rate_of_POC_to_GOC_6
    a₇ = bgc.aggregation_rate_of_POC_to_GOC_7
    a₈ = bgc.aggregation_rate_of_POC_to_GOC_8
    a₉ = bgc.aggregation_rate_of_POC_to_GOC_9

    return sh*a₆*POC^2 + sh*a₇*POC*GOC + a₈*POC*GOC + a₉*POC^2  #39
end

@inline function (bgc::PISCES)(::Val{:POC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust)
    σᶻ = bgc.non_assimilated_fraction.Z
    mᴾ = bgc.phytoplankton_mortality_rate.P
    mᶻ = bgc.zooplankton_quadratic_mortality.Z
    mᴰ = bgc.phytoplankton_mortality_rate.D
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    wₚₒ = bgc.sinking_speed_of_POC
    rᶻ = bgc.zooplankton_linear_mortality.Z
    Kₘ = bgc.half_saturation_const_for_mortality
    b_z, bₘ = bgc.temperature_sensitivity_term
    g_FF = bgc.flux_feeding_rate

    grazing = grazingᶻ(P, D, POC, T, bgc)
    ∑gᶻ = grazing[1]
    gₚₒᶻ = grazing[4]

    sh = get_sh(z, zₘₓₗ)

    R_CaCO₃ = R_CaCO₃(zₘₓₗ, T, P, PAR, bgc) 
    λₚₒ¹ = λ¹(T, O₂, bgc)
    Φ₁ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh, bgc)[1]
    Φ₃ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh, bgc)[3]

    gₚₒᴹ = grazingᴹ(P, D, Z, POC, T, bgc)[4]
    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC #29a
    Φ = Φ(POC, GOC, sh, bgc)

    return σᶻ*∑gᶻ*Z + 0.5*mᴰ*K_mondo(D, Kₘ) + rᶻ*b_z^T*K_mondo(Z, Kₘ)*Z +
     mᶻ*b_z^T*Z^2 + (1 - 0.5*R_CaCO₃)*(mᴾ*K_mondo(P, Kₘ)*P + wᴾ*P^2) + 
     λₚₒ¹*GOC + Φ₁ᴰᴼᶜ + Φ₃ᴰᴼᶜ - (gₚₒᴹ + gₚₒ_FFᴹ)*M - gₚₒᶻ*Z - λₚₒ¹*POC - Φ  #37
end

@inline function (bgc::PISCES)(::Val{:GOC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust)
    σᴹ = bgc.non_assimilated_fraction.M
    mᴾ = bgc.phytoplankton_mortality_rate.P
    mᴰ = bgc.phytoplankton_mortality_rate.D
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    rᴹ = bgc.zooplankton_linear_mortality.M
    Kₘ = bgc.half_saturation_const_for_mortality
    bₘ = bgc.temperature_sensitivity_term.M
    g_FF = bgc.flux_feeding_rate
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms

    ∑gᴹ = grazingᴹ(P, D, Z, POC, T, bgc)[1] 
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    
    Pᵤₚᴹ = Pᵤₚ(M, T, bgc)
    R_CaCO₃ = R_CaCO₃(zₘₓₗ, T, P, PAR, bgc)

    sh = get_sh(z, zₘₓₗ)

    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    wᴰ =  wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ)
    Φ = Φ(POC, GOC, sh, bgc)
    Φ₂ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh, bgc)[2]
    w_GOC = w_GOC(zₑᵤ, zₘₓₗ, bgc)
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC #29b
    λₚₒ¹ = λ¹(T, O₂, bgc)
    

    return σᴹ*(∑gᴹ + ∑g_FFᴹ)*M + rᴹ*bₘ^T*K_mondo(M, Kₘ)*M + Pᵤₚᴹ + 
    0.5*R_CaCO₃*(mᴾ*K_mondo(P, Kₘ)*P + wᴾ*P^2) + 0.5*mᴰ*K_mondo(D, Kₘ)*D^3*wᴰ +
     Φ + Φ₂ᴰᴼᶜ - g_GOC_FFᴹ*M - λₚₒ¹*GOC     #40
end

# This documeent contains functions for:
    # Φ, POC, GOC

@inline function get_Φ(POC, GOC, sh, bgc)
    a₆ = bgc.aggregation_rate_of_POC_to_GOC_6
    a₇ = bgc.aggregation_rate_of_POC_to_GOC_7
    a₈ = bgc.aggregation_rate_of_POC_to_GOC_8
    a₉ = bgc.aggregation_rate_of_POC_to_GOC_9

    return sh*a₆*POC^2 + sh*a₇*POC*GOC + a₈*POC*GOC + a₉*POC^2  #39
end

@inline function (bgc::PISCES)(::Val{:POC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³)
    σᶻ = bgc.non_assimilated_fraction.Z
    mᴾ, mᴰ = bgc.phytoplankton_mortality_rate
    mᶻ = bgc.zooplankton_quadratic_mortality.Z
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    wₚₒ = bgc.sinking_speed_of_POC
    rᶻ = bgc.zooplankton_linear_mortality.Z
    Kₘ = bgc.half_saturation_const_for_mortality
    b_Z, bₘ = bgc.temperature_sensitivity_term
    g_FF = bgc.flux_feeding_rate

    grazing = get_grazingᶻ(P, D, POC, T, bgc)
    ∑gᶻ = grazing[1]
    gₚₒᶻ = grazing[4]

    sh = get_sh(z, zₘₓₗ)

    R_CaCO₃ = get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc) 
    λₚₒ¹ = λ¹(T, O₂, bgc)
    Φ₁ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh, bgc)[1]
    Φ₃ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh, bgc)[3]

    gₚₒᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)[4]
    gₚₒ_FFᴹ = g_FF*(bₘ^T)*wₚₒ*POC #29a
    Φ = get_Φ(POC, GOC, sh, bgc)

    #println("Total change in POC is ", (σᶻ*∑gᶻ*Z + 0.5*mᴰ*concentration_limitation(D, Kₘ)*D + rᶻ*(b_Z^T)*(concentration_limitation(Z, Kₘ) + 3*ΔO₂(O₂, bgc))*Z + mᶻ*(b_Z^T)*Z^2 + (1 - 0.5*R_CaCO₃)*(mᴾ*concentration_limitation(P, Kₘ)*P + sh*wᴾ*P^2) + λₚₒ¹*GOC + Φ₁ᴰᴼᶜ + Φ₃ᴰᴼᶜ - (gₚₒᴹ + gₚₒ_FFᴹ)*M - gₚₒᶻ*Z - λₚₒ¹*POC - Φ))
    #println("Positive terms: a = $(σᶻ*∑gᶻ*Z), b = $(0.5*mᴰ*concentration_limitation(D, Kₘ)*D), c = $(rᶻ*(b_Z^T)*(concentration_limitation(Z, Kₘ) + 3*ΔO₂(O₂, bgc))*Z), " )

    return (σᶻ*∑gᶻ*Z + 0.5*mᴰ*concentration_limitation(D, Kₘ)*D + rᶻ*(b_Z^T)*(concentration_limitation(Z, Kₘ) + 3*ΔO₂(O₂, bgc))*Z 
            + mᶻ*(b_Z^T)*Z^2 + (1 - 0.5*R_CaCO₃)*(mᴾ*concentration_limitation(P, Kₘ)*P + sh*wᴾ*P^2) 
            + λₚₒ¹*GOC + Φ₁ᴰᴼᶜ + Φ₃ᴰᴼᶜ - (gₚₒᴹ + gₚₒ_FFᴹ)*M - gₚₒᶻ*Z - λₚₒ¹*POC - Φ)  #37
end

@inline function (bgc::PISCES)(::Val{:GOC}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³)
    σᴹ = bgc.non_assimilated_fraction.M
    mᴾ = bgc.phytoplankton_mortality_rate.P
    mᴰ = bgc.phytoplankton_mortality_rate.D
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    rᴹ = bgc.zooplankton_linear_mortality.M
    Kₘ = bgc.half_saturation_const_for_mortality
    bₘ = bgc.temperature_sensitivity_term.M
    g_FF = bgc.flux_feeding_rate
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms
    wₚₒ = bgc.sinking_speed_of_POC

    ∑gᴹ = get_grazingᴹ(P, D, Z, POC, T, bgc)[1] 
    ∑g_FFᴹ = get_∑g_FFᴹ(z, zₑᵤ, zₘₓₗ, T, POC, GOC, bgc)
    
    Pᵤₚᴹ = Pᵤₚ(M, T, bgc)
    R_CaCO₃ = get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc)

    sh = get_sh(z, zₘₓₗ)

    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    wᴰ =  wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ)
    Φ = get_Φ(POC, GOC, sh, bgc)
    Φ₂ᴰᴼᶜ = Φᴰᴼᶜ(DOC, POC, GOC, sh, bgc)[2]
    w_GOC = get_w_GOC(z, zₑᵤ, zₘₓₗ, bgc)
    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC #29b
    λₚₒ¹ = λ¹(T, O₂, bgc)

    return (σᴹ*(∑gᴹ + ∑g_FFᴹ)*M + rᴹ*bₘ^T*(concentration_limitation(M, Kₘ) + 3*ΔO₂(O₂, bgc))*M 
            + Pᵤₚᴹ + 0.5*R_CaCO₃*(mᴾ*concentration_limitation(P, Kₘ)*P + sh*wᴾ*P^2) + 0.5*mᴰ*concentration_limitation(D, Kₘ)*D
            + sh*D^2*wᴰ + Φ + Φ₂ᴰᴼᶜ - g_GOC_FFᴹ*M - λₚₒ¹*GOC)     #40 
end
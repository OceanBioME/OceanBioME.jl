#TO DO:
    #How to code ΔCO₃²⁻, as auxiliary field or as in original PISCES?
    #How to define Lₗᵢₘᶜᵃᶜᵒ³()?

#This document contains functions for:
    #R_CaCO₃ (eq77)
    #P_CaCO₃ (eq76)
    #Forcing for CaCO₃ (eq75)

@inline function λ_CaCO₃¹(CaCO₃, bgc, Ω) #no argument required, CaCO₃ given as placeholder
    λ_CaCO₃ = bgc.dissolution_rate_of_calcite
    nca = bgc.exponent_in_the_dissolution_rate_of_calcite
    ΔCO₃²⁻ = max(0, 1 - Ω)
    return λ_CaCO₃*(ΔCO₃²⁻)^nca
end

@inline function get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc) 
    r_CaCO₃ = bgc.rain_ratio_parameter
    Kₙₕ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_ammonium.P
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton.P
    Pₘₐₓ = bgc.threshold_concentration_for_size_dependency.P
    P₁ =  I₁(P, Pₘₐₓ)
    P₂ = I₂(P, Pₘₐₓ)
    Lₙᴾ = P_nutrient_limitation(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, bgc)[5]
    Kₙₕ₄ᴾ = Kᵢᴶ(Kₙₕ₄ᴾᵐⁱⁿ, P₁, P₂, Sᵣₐₜᴾ)
    Lₗᵢₘᶜᵃᶜᵒ³ = min(Lₙᴾ, concentration_limitation(Fe, 6e-11), concentration_limitation(PO₄, Kₙₕ₄ᴾ))
    return r_CaCO₃*Lₗᵢₘᶜᵃᶜᵒ³*T*max(1, P/2)*max(0, PAR - 1)*30*(1 + exp((-(T-10)^2)/25))*min(1, 50/(zₘₓₗ + eps(0.0)))/((0.1 + T)*(4 + PAR)*(30 + PAR)) #eq77
end

@inline function P_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, D, Z, M, POC, T, PAR, zₘₓₗ, z, bgc) 
    mᴾ = bgc.phytoplankton_mortality_rate.P
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    ηᶻ = bgc.proportion_of_sinking_grazed_shells.Z
    ηᴹ = bgc.proportion_of_sinking_grazed_shells.M
    sh = get_sh(z, zₘₓₗ)
    
    return get_R_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, T, PAR, zₘₓₗ, bgc)*(ηᶻ*get_grazingᶻ(P, D, POC, T, bgc)[2]*Z+ηᴹ*get_grazingᴹ(P, D, Z, POC, T, bgc)[2]*M + 0.5*(mᴾ*concentration_limitation(P, Kₘ)*P + sh*wᴾ*P^2)) #eq76
end

@inline function (bgc::PISCES)(::Val{:CaCO₃}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR¹, PAR², PAR³)

    return P_CaCO₃(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, Fe, D, Z, M, POC, T, PAR, zₘₓₗ, z, bgc) - λ_CaCO₃¹(CaCO₃, bgc, Ω)*CaCO₃ #partial derivative omitted as sinking is accounted for in other parts of model
end
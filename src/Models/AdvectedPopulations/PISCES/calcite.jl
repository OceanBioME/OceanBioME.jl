#TO DO:
    #How to code ΔCO₃²⁻, as auxiliary field or as in original PISCES?
    #How to define Lₗᵢₘᶜᵃᶜᵒ³()?

#This document contains functions for:
    #R_CaCO₃ (eq77)
    #P_CaCO₃ (eq76)
    #Forcing for CaCO₃ (eq75)

@inline function λ_CaCO₃¹(CaCO₃, bgc) #no argument required, CaCO₃ given as placeholder
    λ_CaCO₃ = bgc.dissolution_rate_of_calcite
    nca = bgc.exponent_in_the_dissolution_rate_of_calcite
    #Ω = bgc.carbonate_sat_ratio #define this as an auxiliary field, or using Nemo source code as in PISCES?
    Ω = 0
    ΔCO₃²⁻ = max(0, 1 - Ω)
    return λ_CaCO³*(ΔCO₃²⁻)^nca
end

@inline function R_CaCO₃(P, T, PAR, zₘₓₗ, bgc) 
    r_CaCO₃ = bgc.rain_ratio_parameter
    Lₗᵢₘᶜᵃᶜᵒ³ = #does this equal 1 or as defined in original PISCES?
    return r_CaCO₃*Lₗᵢₘᶜᵃᶜᵒ³*T*max(1, P/2)*max(0, PAR - 1)*30*(1 + exp((-(T-10)^2)/25))*min(1, 50/zₘₓₗ)/((0.1 + T)*(4 + PAR)*(30 + PAR)) #eq77
end

@inline function P_CaCO₃(P, Z, M, T, PAR, zₘₓₗ, z, bgc) 
    mᴾ = bgc.zooplankton_quadratic_mortality.P
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    ηᶻ = bgc.proportion_of_sinking_grazed_shells.Z
    ηᴹ = bgc.proportion_of_sinking_grazed_shells.M
    sh = get_sh(z, zₘₓₗ)
    
    return R_CaCO₃(P, T, PAR, zₘₓₗ, bgc)*(ηᶻ*grazingᶻ(P, D, POC, T, bgc)[2]*Z+ηᴹ*grazingᴹ(P, D, Z, POC, T, bgc)[2]*M + 0.5*(mᴾ*K_mondo(P, Kₘ)*P + sh*wᴾ*P^2)) #eq76
end

@inline function (bgc::PISCES)(::Val{:CaCO₃}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust) 

    return P_CaCO₃(P, Z, M, T, PAR, zₘₓₗ, z, bgc) - λ_CaCO₃¹(CaCO₃, bgc)*CaCO₃ #partial derivative omitted as sinking is accounted for in other parts of model
end
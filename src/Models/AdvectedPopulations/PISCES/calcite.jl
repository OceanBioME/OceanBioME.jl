#TO DO:
    #Write a partial derivative in (75)
    #How to code ΔCO₃²⁻, where do we get the value of CO₃²⁻ from?
    #What is η? - check original PISCES documentation.
    #What is Lₗᵢₘᶜᵃᶜᵒ³()?

#This document contains functions for:
    #R_CaCO₃ (eq77)
    #P_CaCO₃ (eq76)
    #Forcing for CaCO₃ (eq75)

@inline function λ_CaCO₃¹()
    λ_CaCO₃ = bgc.dissolution_rate_of_calcite
    nca = bgc.exponent_in_the_dissolution_rate_of_calcite
    ΔCO₃²⁻ = 0 #how to define this?
    return λ_CaCO³*(ΔCO₃²⁻)^nca
end

@inline function R_CaCO₃(zₘₓₗ, T, P, PAR) 
    r_CaCO₃ = bgc.rain_ratio_parameter
    return r_CaCO₃*Lₗᵢₘᶜᵃᶜᵒ³()*T*max(1, P/2)*max(0, PAR - 1)*30*(1 + exp((-(T-10)^2)/25))*min(1, 50/zₘₓₗ)/((0.1 + T)*(4 + PAR)*(30 + PAR)) #eq77
end

@inline function P_CaCO₃(zₘₓₗ, T, PAR, P, M) 
    mᴾ = bgc.zooplankton_quadratic_mortality[1]
    Kₘ = bgc.half_saturation_const_for_mortality
    ωᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    return R_CaCO₃()*(ηᶻ*grazingᶻ()[2]*Z+ηᴹ*grazingᴹ[2]*M + 0.5*(mᴾ*K_mondo(P, Kₘ)*P) + sh*ωᴾ*P^2) #eq76
end

@inline function λ_CaCO₃¹()
    ΔCO₃²⁻ = #Ask Jago
    return 0
end

@inline function (pisces::PISCES)(::Val{:CaCO₃}, x, y, z, t, P, T, PAR) 

    return P_CaCO₃() - λ_CaCO₃¹()*CaCO₃ - ω_goc* #(75) how to write partial derivative here
end
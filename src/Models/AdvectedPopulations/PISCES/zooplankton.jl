# Grazing function is defined below. Similar to the L parameter, it returns all the components as they get used frequently
# Jₜₕᵣₑₛₕᶻ seems to only be defined for micro and meso zooplankton rather than each species.

@inline function grazingᶻ(P, D, POC, T) 
    pₚᶻ = bgc.preference_for_nanophytoplankton[1]
    p_Dᶻ = bgc.preference_for_diatoms[1]
    pₚₒᶻ = bgc.preference_for_POC[1]
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Fₜₕᵣₑₛₕᶻ = bgc.food_threshold_for_zooplankton[1]
    gₘₐₓᶻ = bgc.max_grazing_rate[1]
    K_Gᶻ = bgc.half_saturation_const_for_grazing[1]

    F = pₚᶻ*max(0, P - Jₜₕᵣₑₛₕᶻ) + p_Dᶻ*max(0, D - Jₜₕᵣₑₛₕᶻ) + pₚₒᶻ*max(0, POC - Jₜₕᵣₑₛₕᶻ)
    Fₗᵢₘ = max(0, F - min(0.5*F, Fₜₕᵣₑₛₕᶻ))

    grazing_arg = gₘₐₓᶻ*b_z^T*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᶻ + pₚᶻ*P + p_Dᶻ*D + pₚₒᶻ*POC + eps(0.0)))

    gₚᶻ = (pₚᶻ*max(0, P - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    g_Dᶻ = (p_Dᶻ*max(0, D - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    gₚₒᶻ = (pₚₒᶻ*max(0, POC - Jₜₕᵣₑₛₕᶻ))*grazing_arg #26a
    
    ∑gᶻ= gₚᶻ + g_Dᶻ + gₚₒᶻ  #Sum grazing rates on each prey species for microzooplankton

    return ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ
end

@inline function grazingᴹ(P, D, Z, POC, T) 
    pₚᴹ = bgc.preference_for_nanophytoplankton[2]
    p_Dᴹ = bgc.preference_for_diatoms[2]
    pₚₒᴹ = bgc.preference_for_POC[2]
    p_zᴹ = bgc.preference_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    Fₜₕᵣₑₛₕᴹ = bgc.food_threshold_for_zooplankton[2]
    gₘₐₓᴹ = bgc.max_grazing_rate[2]
    K_Gᴹ = bgc.half_saturation_const_for_grazing[2]
    
    F = pₚᴹ*max(0, P - Jₜₕᵣₑₛₕᴹ) + p_Dᴹ*max(0, D - Jₜₕᵣₑₛₕᴹ) + pₚₒᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ) + p_zᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ)
    Fₗᵢₘ = max(0, F - min(0.5*F, Fₜₕᵣₑₛₕᴹ))

    grazing_arg =  gₘₐₓᴹ*bₘ^T*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᴹ + pₚᴹ*P + p_Dᴹ*D + pₚₒᴹ*POC + p_zᴹ*Z + eps(0.0)))

    gₚᴹ = (pₚᴹ*max(0, P - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    g_Dᴹ = (p_Dᴹ*max(0, D - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    gₚₒᴹ = (pₚₒᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    g_zᴹ = (p_zᴹ*max(0, Z - Jₜₕᵣₑₛₕᴹ))*grazing_arg #26a
    ∑gᴹ = gₚᴹ +  g_Dᴹ + gₚₒᴹ + g_Zᴹ #Sum grazing rates on each prey species for mesozooplankton
    
    return  ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ 
end

@inline function w_GOC(zₑᵤ, zₘₓₗ)
    zₘₐₓ = max(zₑᵤ, zₘₓₗ) 
    w_GOCᵐⁱⁿ = bgc.min_sinking_speed_of_GOC
    return w_GOCᵐⁱⁿ + (200 - w_GOCᵐⁱⁿ)*(max(0, z-zₘₐₓ))/(5000) #41b
end

@inline function ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)
    wₚₒ = bgc.sinking_speed_of_POC
    g_FF = bgc.flux_feeding_rate
    bₘ = bgc.temperature_sensitivity_term[2]

    w_GOC = w_GOC(zₑᵤ, zₘₓₗ)

    gₚₒ_FFᴹ = g_FF*bₘ^T*wₚₒ*POC #29a
    g_GOC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC #29b
    return g_GOC_FFᴹ + gₚₒ_FFᴹ
end

# gross growth efficiency, defined for both but g_zᴹ and Z do not appear for eᶻ so have passed in as 0 and 1 respectively to avoid divide by zero error.
@inline function eᴶ(eₘₐₓᴶ, σᴶ, gₚᴶ, g_Dᴶ, gₚₒᴶ, g_zᴹ, N, Fe, P, D, POC, Z, J)
    θᴺᶜ = bgc.NC_redfield_ratio

    ∑ᵢθᴺᴵgᵢᴶ = θ(N,P)*gₚᴶ + θ(N, D)*g_Dᴶ + θ(N, POC)*gₚₒᴶ + θ(N, Z)*g_zᴹ
    ∑ᵢθᶠᵉᴵgᵢᴶ = θ(Fe, P)*gₚᴶ + θ(Fe, D)*g_Dᴶ + θ(Fe, POC)*gₚₒᴶ + θ(Fe, Z)*g_zᴹ
    ∑ᵢgᵢᴶ = gₚᴶ + g_Dᴶ + gₚₒᴶ + g_zᴹ

    eₙᴶ = min(1, (∑ᵢθᴺᴵgᵢᴶ)/(θᴺᶜ*∑ᵢgᵢᴶ), (∑ᵢθᶠᵉᴵgᵢᴶ)/(θ(Fe, Z)*∑ᵢgᵢᴶ))   #27a

    return eₙᴶ*min(eₘₐₓᴶ, (1 - σᴶ)* (∑ᵢθᶠᵉᴵgᵢᴶ)/(θ(Fe, J)*∑ᵢgᵢᴶ)) #27b
end


@inline function (pisces::PISCES)(::Val{:Z}, x, y, z, t, Z, M, P, POC, D, T, O₂)    #args not correct
    mᶻ = bgc.zooplankton_quadratic_mortality[1]
    b_z = bgc.temperature_sensitivity_term
    Kₘ = bgc.half_saturation_const_for_mortality
    rᶻ = bgc.zooplankton_linear_mortality[1]
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton[1]
    σᶻ = bgc.non_assimilated_fraction[1]

    ∑gᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ = grazingᶻ(P, D, POC, T) 
    g_zᴹ = grazingᴹ(P, D, Z, POC, T)[5]

    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, N, Fe, P, D, POC, 1, Z)

    return eᶻ*(gₚᶻ + g_Dᶻ + gₚₒᶻ)*Z - g_zᴹ*M - mᶻ*b_z^T*Z^2 - rᶻ*b_z^T*(K_mondo(Z, Kₘ) + 3*ΔO₂(O₂))*Z   #24
end

@inline function (pisces::PISCES)(::Val{:M}, x, y, z, t, Z, M, P, POC, GOC, D, T, O₂, zₘₓₗ, zₑᵤ) #args not correct
    mᴹ = bgc.zooplankton_quadratic_mortality[2]
    bₘ = bgc.temperature_sensitivity_term[2]
    rᴹ = bgc.zooplankton_linear_mortality[2]
    Kₘ = bgc.half_saturation_const_for_mortality

    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton[2]
    σᴹ = bgc.non_assimilated_fraction[2]

    ∑gᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ  = grazingᴹ(P, D, Z, POC, T) 

    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)
    
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ, N, Fe, P, D, POC, Z, M)

    return eᴹ*(gₚᴹ + g_Dᴹ + gₚₒᴹ + ∑g_FFᴹ + g_zᴹ)*M - mᴹ*bₘ^T*M^2 - rᴹ*bₘ^T*(K_mondo(M, Kₘ) + 3*ΔO₂(O₂))*M   #28
end
@inline function grazing_argᶻ(P, POC, D, T)
    pₚᶻ = bgc.preference_for_nanophytoplankton[1]
    p_Dᶻ = bgc.preference_for_diatoms[1]
    pₚₒᶻ = bgc.preference_for_POC[1]
    b_z = bgc.temperature_sensitivity_term[1]
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton    #Think this can be changed for each species
    Fₜₕᵣₑₛₕᶻ = bgc.food_threshold_for_zooplankton[1]
    gₘₐₓᶻ = bgc.max_grazing_rate[1]
    K_Gᶻ = bgc.half_saturation_const_for_grazing[1]
    
    F = pₚᶻ*max(0, P - Jₜₕᵣₑₛₕᶻ) + p_Dᶻ*max(0, D - Jₜₕᵣₑₛₕᶻ) + pₚₒᶻ*max(0, POC - Jₜₕᵣₑₛₕᶻ)
    Fₗᵢₘ = max(0, F - min(0.5*F, Fₜₕᵣₑₛₕᶻ))
    
    return gₘₐₓᶻ*b_z^T*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᶻ + pₚᶻ*P + p_Dᶻ*D + pₚₒᶻ*POC + eps(0.0)))
end

@inline function grazing_argᴹ(P, POC, D, Z, T)
    pₚᴹ = bgc.preference_for_nanophytoplankton[2]
    p_Dᴹ = bgc.preference_for_diatoms[2]
    pₚₒᴹ = bgc.preference_for_POC[2]
    p_zᴹ = bgc.preference_for_microzooplankton
    bₘ = bgc.temperature_sensitivity_term[2]
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton   #Think this can be changed for each species
    Fₜₕᵣₑₛₕᴹ = bgc.food_threshold_for_zooplankton[2]
    gₘₐₓᴹ = bgc.max_grazing_rate[2]
    K_Gᴹ = bgc.half_saturation_const_for_grazing[2]

    F = pₚᴹ*max(0, P - Jₜₕᵣₑₛₕᴹ) + p_Dᴹ*max(0, D - Jₜₕᵣₑₛₕᴹ) + pₚₒᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ) + p_zᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ)
    Fₗᵢₘ = max(0, F - min(0.5*F, Fₜₕᵣₑₛₕᴹ))

    return gₘₐₓᴹ*bₘ^T*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᴹ + pₚᴹ*P + p_Dᴹ*D + pₚₒᴹ*POC + p_zᴹ*Z + eps(0.0)))

end


@inline gᴶ(I, pᵢᴶ, Iₜₕᵣₑₛₕᴶ, grazing_arg) = (pᵢᴶ*max(0, I - Iₜₕᵣₑₛₕᴶ))*grazing_arg #26a

@inline function ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)
    w_POC = bgc.sinking_speed_of_POC
    g_FF = bgc.flux_feeding_rate
    w_GOCᵐⁱⁿ = bgc.min_sinking_speed_of_GOC
    bₘ = bgc.temperature_sensitivity_term[2]

    zₘₐₓ = max(zₑᵤ, zₘₓₗ)   #41a
    w_GOC = w_GOCᵐⁱⁿ + (200 - w_GOCᵐⁱⁿ)*(max(0, z-zₘₐₓ))/(5000) #41b

    g_GOC_FFᴹ = g_FF*bₘ^T*w_POC*POC #29b
    g_POC_FFᴹ = g_FF*bₘ^T*w_GOC*GOC #29a
    return g_GOC_FFᴹ + g_POC_FFᴹ
end

@inline function eᴶ(eₘₐₓᴶ, σᴶ, gₚᴶ, g_Dᴶ, gₚₒᴶ, g_zᴹ, N, Fe, P, D, POC, Z, J)
    θᴺᶜ = bgc.NC_redfield_ratio

    ∑ᵢθᴺᴵgᵢᴶ = θ(N,P)*gₚᴶ + θ(N, D)*g_Dᴶ + θ(N, POC)*gₚₒᴶ + θ(N, Z)*g_zᴹ
    ∑ᵢθᶠᵉᴵgᵢᴶ = θ(Fe, P)*gₚᴶ + θ(Fe, D)*g_Dᴶ + θ(Fe, POC)*gₚₒᴶ + θ(Fe, Z)*g_zᴹ
    ∑ᵢgᵢᴶ = gₚᴶ + g_Dᴶ + gₚₒᴶ + g_zᴹ

    eₙᴶ = min(1, (∑ᵢθᴺᴵgᵢᴶ)/(θᴺᶜ*∑ᵢgᵢᴶ), (∑ᵢθᶠᵉᴵgᵢᴶ)/(θ(Fe, Z)*∑ᵢgᵢᴶ))   #27a

    return eₙᴶ*min(eₘₐₓᴶ, (1 - σᴶ)* (∑ᵢθᶠᵉᴵgᵢᴶ)/(θ(Fe, J)*∑ᵢgᵢᴶ)) #27b
end


@inline function (pisces::PISCES)(::Val{:Z}, x, y, z, t, Z, M, P, POC, D, T, O₂) 
    mᶻ = bgc.zooplankton_quadratic_mortality[1]
    b_z = bgc.temperature_sensitivity_term
    Kₘ = bgc.half_saturation_const_for_mortality
    rᶻ = bgc.zooplankton_linear_mortality[1]
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton[1]
    σᶻ = bgc.non_assimilated_fraction[1]
    pₚᶻ = bgc.preference_for_nanophytoplankton[1]
    p_Dᶻ = bgc.preference_for_diatoms[1]
    pₚₒᶻ = bgc.preference_for_POC[1]
    p_zᴹ = bgc.preference_for_microzooplankton
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton

    
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)

    gₚᶻ = gᴶ(P, pₚᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    g_Dᶻ = gᴶ(D, p_Dᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    gₚₒᶻ = gᴶ(POC, pₚₒᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    g_zᴹ = gᴶ(Z, p_zᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)

    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, 0, N, Fe, P, D, POC, 1, Z)

    return eᶻ*(gₚᶻ + g_Dᶻ + gₚₒᶻ)*Z - g_zᴹ*M - mᶻ*b_z^T*Z^2 - rᶻ*b_z^T*(K_mondo(Z, Kₘ) + 3*ΔO₂(O₂))*Z   #24
end

@inline function (pisces::PISCES)(::Val{:M}, x, y, z, t, Z, M, P, POC, GOC, D, T, O₂, zₘₓₗ, zₑᵤ)
    mᴹ = bgc.zooplankton_quadratic_mortality[2]
    bₘ = bgc.temperature_sensitivity_term[2]
    rᴹ = bgc.zooplankton_linear_mortality[2]
    Kₘ = bgc.half_saturation_const_for_mortality

    eₘₐₓᴹ = bgc.max_growth_efficiency_of_zooplankton[2]
    σᴹ = bgc.non_assimilated_fraction[2]

    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    pₚᴹ = bgc.preference_for_nanophytoplankton[2]
    p_Dᴹ = bgc.preference_for_diatoms[2]
    pₚₒᴹ = bgc.preference_for_POC[2]
    p_zᴹ = bgc.preference_for_microzooplankton

    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)

    gₚᴹ = gᴶ(P, pₚᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)
    g_Dᴹ = gᴶ(D, p_Dᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)
    gₚₒᴹ = gᴶ(POC, pₚₒᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)
    g_zᴹ = gᴶ(Z, p_zᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)

    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)
    
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ, N, Fe, P, D, POC, Z, M)

    return eᴹ*(gₚᴹ + g_Dᴹ + gₚₒᴹ + ∑g_FFᴹ + g_zᴹ)*M - mᴹ*bₘ^T*M^2 - rᴹ*bₘ^T*(K_mondo(M, Kₘ) + 3*ΔO₂(O₂))*M   #28
end
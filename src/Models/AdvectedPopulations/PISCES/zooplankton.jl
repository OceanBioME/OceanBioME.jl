@inline function grazing_argᶻ(P, POC, D, T)
    pₚᶻ = bgc.preference_for_nanophytoplankton[1]
    p_Dᶻ = bgc.preference_for_diatoms[1]
    pₚₒᶻ = bgc.preference_for_POC[1]
    b_z = bgc.temperature_sensitivity_term
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
    b_z = bgc.temperature_sensitivity_term
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton   #Think this can be changed for each species
    Fₜₕᵣₑₛₕᴹ = bgc.food_threshold_for_zooplankton[2]
    gₘₐₓᴹ = bgc.max_grazing_rate[2]
    K_Gᴹ = bgc.half_saturation_const_for_grazing[2]

    F = pₚᴹ*max(0, P - Jₜₕᵣₑₛₕᴹ) + p_Dᴹ*max(0, D - Jₜₕᵣₑₛₕᴹ) + pₚₒᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ) + p_zᴹ*max(0, POC - Jₜₕᵣₑₛₕᴹ)
    Fₗᵢₘ = max(0, F - min(0.5*F, Fₜₕᵣₑₛₕᴹ))

    return gₘₐₓᴹ*b_z^T*(Fₗᵢₘ)/((F + eps(0.0))*(K_Gᴹ + pₚᴹ*P + p_Dᴹ*D + pₚₒᴹ*POC + p_zᴹ*Z + eps(0.0)))

end


@inline gᴶ(I, pᵢᴶ, Iₜₕᵣₑₛₕᴶ, grazing_arg) = (pᵢᴶ*max(0, I - Iₜₕᵣₑₛₕᴶ))*grazing_arg #26a



@inline function eᴶ(eₘₐₓᴶ, σᴶ, gₚᴶ, g_Dᴶ, gₚₒᴶ, N, P, D, POC, J)
    θᴺᶜ = bgc.NC_redfield_ratio

    Σᵢθᴺᴵgᵢᴶ = θ(N,P)*gₚᴶ + θ(N, D)*g_Dᴶ + θ(N, POC)*gₚₒᴶ 
    Σᵢθᶠᵉᴵgᵢᴶ = θ(Fe, P)*gₚᴶ + θ(Fe, D)*g_Dᴶ + θ(Fe, POC)*gₚₒᴶ 
    Σᵢgᵢᴶ = gₚᴶ + g_Dᴶ + gₚₒᴶ 

    eₙᴶ = min(1, (Σᵢθᴺᴵgᵢᴶ)/(θᴺᶜ*Σᵢgᵢᴶ), (Σᵢθᶠᵉᴵgᵢᴶ)/(θ(Fe, Z)*Σᵢgᵢᴶ))   #27a

    return eₙᴶ*min(eₘₐₓᴶ, (1 - σᴶ)* (Σᵢθᶠᵉᴵgᵢᴶ)/(θ(Fe, J)*Σᵢgᵢᴶ)) #27b
end


@inline function (pisces::PISCES)(::Val{:Z}, x, y, z, t, Z, T, O₂) 
    mᶻ = bgc.zooplankton_quadratic_mortality[1]
    b_z = bgc.temperature_sensitivity_term
    Kₘ = bgc.half_saturation_const_for_mortality
    rᶻ = bgc.zooplankton_linear_mortality[1]
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton[1]
    σᶻ = bgc.non_assimilated_fraction[1]
    pₚᶻ = bgc.preference_for_nanophytoplankton[1]
    p_Dᶻ = bgc.preference_for_diatoms[1]
    pₚₒᶻ = bgc.preference_for_POC[1]
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton

    
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)

    gₚᶻ = gᴶ(P, pₚᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    g_Dᶻ = gᴶ(D, p_Dᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    gₚₒᶻ = gᴶ(POC, pₚₒᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    g_zᴹ = gᴶ(, pₚᶻ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)

    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, N, P, D, POC, Z)

    return eᶻ*(gₚᶻ + g_Dᶻ + gₚₒᶻ)*Z - g_zᴹ*M - mᶻ*b_z^T*Z^2 - rᶻ*b_z^T*(K_mondo(Z, Kₘ) + 3*ΔO₂(O₂))*Z   #24
end

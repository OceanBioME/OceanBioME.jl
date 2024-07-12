@inline function eᴵ(eₘₐₓᴵ, σᴵ)
    θᴺᶜ = bgc.NC_redfield_ratio

    Σᵢθᴺᴵgᵢᶻ = 
    Σᵢθᶠᵉᴵgᵢᶻ = 
    Σᵢgᵢᶻ = 

    eₙᴵ = min(1, (Σᵢθᴺᴵgᵢᶻ)/(θᴺᶜ*Σᵢgᵢᶻ), (Σᵢθᶠᵉᴵgᵢᶻ)/(θ(Fe, Z)*Σᵢgᵢᶻ))   #27a

    return eₙᴵ*min(eₘₐₓᴵ, (1 - σᴵ)* (Σᵢθᶠᵉᴵgᵢᶻ)/(θ(Fe, Z)*Σᵢgᵢᶻ)) #27b
end


@inline function (pisces::PISCES)(::Val{:Z}, x, y, z, t, Z, T, O₂) 
    mᶻ = bgc.zooplankton_quadratic_mortality[1]
    b_z = bgc.temperature_sensitivity_term
    Kₘ = bgc.half_saturation_const_for_mortality
    rᶻ = bgc.zooplankton_linear_mortality[1]
    eₘₐₓᶻ = bgc.max_growth_efficiency_of_zooplankton[1]
    σᶻ = bgc.non_assimilated_fraction[1]

    eᶻ = eᴵ(eₘₐₓᶻ, σᶻ)
    gₚᶻ = 
    g_Dᶻ = 
    gₚₒᶻ = 
    g_zᴹ = 

    return eᶻ*(gₚᶻ + g_Dᶻ + gₚₒᶻ)*Z - g_zᴹ*M - mᶻ*b_z^T*Z^2 - rᶻ*b_z^T*(K_mondo(Z, Kₘ) + 3*ΔO₂(O₂))*Z   #24
end

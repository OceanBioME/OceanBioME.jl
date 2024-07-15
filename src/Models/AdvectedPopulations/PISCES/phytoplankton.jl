@inline θ(I,J) = I/(J + eps(0.0))   #eq 0
@inline K_mondo(I, J) = I/(I + J + eps(0.0))
@inline Cₚᵣₒ(I, Iᶜʰˡ, PARᴵ, L_day, αᴵ, μₚ, Lₗᵢₘᴵ)=1-exp(-αᴵ*(θ(Iᶜʰˡ,I))*PARᴵ/(L_day*μₚ*Lₗᵢₘᴵ + eps(0.0)))

@inline f₁(L_day) = 1.5*L_day/(0.5+L_day)  #eq 3a
@inline t_dark(zₘₓₗ, zₑᵤ) = max(0, zₘₓₗ-zₑᵤ)^2/86400#eq 3b,c
@inline f₂(zₘₓₗ, zₑᵤ, t_darkᴵ) = 1 - t_dark(zₘₓₗ, zₑᵤ)/(t_dark(zₘₓₗ, zₑᵤ) + t_darkᴵ) #eq 3d

@inline fₚ(T) = bgc.temperature_sensitivity_of_growth^T #eq 4a

@inline L_NH₄(NO₃, NH₄, Kₙₒ₃ᴵ, Kₙₕ₄ᴵ) = Kₙₒ₃ᴵ*NH₄/(Kₙₒ₃ᴵ*Kₙₕ₄ᴵ+Kₙₕ₄ᴵ*NO₃+Kₙₒ₃ᴵ*NH₄) #eq 6d
@inline L_NO₃(NO₃, NH₄, Kₙₒ₃ᴵ, Kₙₕ₄ᴵ) = Kₙₕ₄ᴵ*NO₃/(Kₙₒ₃ᴵ*Kₙₕ₄ᴵ+Kₙₕ₄ᴵ*NO₃+Kₙₒ₃ᴵ*NH₄) #eq 6e
@inline L_Fe(I, Iᶠᵉ, θₒₚₜᶠᵉᴵ, θₘᵢₙᶠᵉᴵ) = min(1, max(0, (θ(Iᶠᵉ, I) - θₘᵢₙᶠᵉᴵ)/θₒₚₜᶠᵉᴵ)) #eq 6f

@inline θᶠᵉₘᵢₙ(I, Iᶜʰˡ, Lₙᴵ, Lₙₒ₃ᴵ) = 0.0016/(55.85) * θ(Iᶜʰˡ, I) + 1.21e-5*14*Lₙᴵ/(55.85*7.625)*1.5+1.15*14*Lₙₒ₃ᴵ/(55.85*7.625) #eq 20 -> Lₙ could be meant to be L_NH₄?

@inline I₁(I, Iₘₐₓ) = min(I, Iₘₐₓ) #eq 7a
@inline I₂(I, Iₘₐₓ) = max(0, I - Iₘₐₓ) #eq 7b
@inline Kᵢᴶ(Kᵢᴶᵐⁱⁿ, J₁, J₂, Sᵣₐₜᴶ) = Kᵢᴶᵐⁱⁿ* (J₁ + Sᵣₐₜᴶ* J₂)/(J₁ + J₂) #eq 7c

@inline function μᴵᶠᵉ(I, Iᶠᵉ, θₘₐₓᶠᵉᴵ, Sᵣₐₜᴵ, K_Feᴵᶠᵉᵐⁱⁿ, Iₘₐₓ, L_Feᴵ, bFe)
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero

    μₚ = μ⁰ₘₐₓ*fₚ(T) #4b

    I₂ = max(0, I - Iₘₐₓ) #18c
    I₁ = I - I₂     #18c

    K_Feᴵᶠᵉ = K_Feᴵᶠᵉᵐⁱⁿ*(I₁ + Sᵣₐₜᴵ*I₂)/(I₁+I₂)    #18b

    Lₗᵢₘ₁ᴵᶠᵉ = K_mondo(bFe, K_Feᴵᶠᵉ)    #18a
    Lₗᵢₘ₂ᴵᶠᵉ = (4 - 4.5*L_Feᴵ)/(L_Feᴵ + 0.5) #19


    return θₘₐₓᶠᵉᴵ*Lₗᵢₘ₁ᴵᶠᵉ*Lₗᵢₘ₂ᴵᶠᵉ*(1 - (θ(Iᶠᵉ, I))/(θₘₐₓᶠᵉᴵ))/(1.05 - (θ(Iᶠᵉ, I))/(θₘₐₓᶠᵉᴵ))*μₚ  #17
end

@inline function μᴵ(I, Iᶜʰˡ, PARᴵ, L_day, T, αᴵ, Lₗᵢₘᴵ, zₘₓₗ, zₑᵤ, t_darkᴵ)
    
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero

    μₚ = μ⁰ₘₐₓ*fₚ(T)    #eq 4b      #address t_darkᴵ

    return μₚ * f₁(L_day) * f₂(zₘₓₗ, zₑᵤ, t_darkᴵ) * Cₚᵣₒ(I, Iᶜʰˡ, PARᴵ, L_day, αᴵ, μₚ, Lₗᵢₘᴵ) * Lₗᵢₘᴵ #2b 
end

@inline function Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)
    θₒₚₜᶠᵉᵖ = bgc.optimal_iron_quota[1]
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton[1]
    Kₙₒ₃ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_nitrate[1]
    Kₙₕ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_ammonium[1]
    Pₘₐₓ =

    P₁ = I₁(P, Pₘₐₓ)
    P₂ = I₂(P, Pₘₐₓ)

    Kₙₒ₃ᴾ = Kᵢᴶ(Kₙₒ₃ᴾᵐⁱⁿ, P₁, P₂, Sᵣₐₜᴾ)
    Kₙₕ₄ᴾ = Kᵢᴶ(Kₙₕ₄ᴾᵐⁱⁿ, P₁, P₂, Sᵣₐₜᴾ)

    Lₚₒ₄ᴾ = K_mondo(PO₄, Kₚₒ₄ᴾ) #6b
    Lₙₕ₄ᴾ = L_NH₄(NO₃, NH₄, Kₙₒ₃ᴾ, Kₙₕ₄ᴾ)
    Lₙₒ₃ᴾ = L_NO₃(NO₃, NH₄, Kₙₒ₃ᴾ, Kₙₕ₄ᴾ)
    Lₙᴾ = Lₙₒ₃ᴾ + Lₙₕ₄ᴾ         #6c

    θₘᵢₙᶠᵉᵖ = θᶠᵉₘᵢₙ(P, Pᶜʰˡ, Lₙᴾ, Lₙₒ₃ᴾ)
    L_Feᴾ = L_Fe(P, Pᶠᵉ, θₒₚₜᶠᵉᵖ, θₘᵢₙᶠᵉᵖ)

    return min(Lₚₒ₄ᴾ, Lₙᴾ, L_Feᴾ), Lₚₒ₄ᴾ, Lₙₕ₄ᴾ, Lₙₒ₃ᴾ, Lₙᴾ, L_Feᴾ #6a
end

@inline function Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)
    θₒₚₜᶠᵉᴰ = bgc.optimal_iron_quota[2]
    Sᵣₐₜᴰ = bgc.size_ratio_of_phytoplankton[2]
    Kₙₒ₃ᴰᵐⁱⁿ = bgc.min_half_saturation_const_for_nitrate[2]
    Kₙₕ₄ᴰᵐⁱⁿ = bgc.min_half_saturation_const_for_ammonium[2]
    Kₛᵢᴰᵐⁱⁿ = bgc.min_half_saturation_const_for_silicate
    Kₛᵢ = bgc.parameter_for_half_saturation_const
    Dₘₐₓ = 
    SI = #Si with a crescent moon on it

    D₁ = I₁(D, Dₘₐₓ)
    D₂ = I₂(D, Dₘₐₓ)

    Kₙₒ₃ᴰ = Kᵢᴶ(Kₙₒ₃ᴰᵐⁱⁿ, D₁, D₂, Sᵣₐₜᴰ)
    Kₙₕ₄ᴰ = Kᵢᴶ(Kₙₕ₄ᴰᵐⁱⁿ, D₁, D₂, Sᵣₐₜᴰ)

    Lₚₒ₄ᴰ = K_mondo(PO₄, Kₚₒ₄ᴰ) #6b
    Lₙₕ₄ᴰ = L_NH₄(NO₃, NH₄, Kₙₒ₃ᴰ, Kₙₕ₄ᴰ)
    Lₙₒ₃ᴰ = L_NO₃(NO₃, NH₄, Kₙₒ₃ᴰ, Kₙₕ₄ᴰ)
    Lₙᴰ = Lₙₒ₃ᴰ + Lₙₕ₄ᴰ         #6c

    θₘᵢₙᶠᵉᴰ = θᶠᵉₘᵢₙ(D, Dᶜʰˡ, Lₙᴰ, Lₙₒ₃ᴰ)
    L_Feᴰ = L_Fe(D, Dᶠᵉ ,θₒₚₜᶠᵉᴰ, θₘᵢₙᶠᵉᴰ)
    Kₛᵢᴰ = Kₛᵢᴰᵐⁱⁿ + 7*SI^2 / (Kₛᵢ^2 + SI^2) #12
    Lₛᵢᴰ = K_mondo(Si, Kₛᵢᴰ)    #11b

    return min(Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ, Lₛᵢᴰ), Lₚₒ₄ᴰ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ, Lₙᴰ, Lₛᵢᴰ, L_Feᴰ    #11a
end

@inline function fθₒₚₜˢⁱᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, μᴰ, T, ϕ)
    θₘˢⁱᴰ = bgc.optimal_SiC_uptake_ratio_of_diatoms
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero
    Kₛᵢ¹ = bgc.parameter_for_SiC[1]
    Kₛᵢ² = bgc.parameter_for_SiC[2]

    Lₗᵢₘᴰ, Lₚₒ₄ᴰ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ, Lₙᴰ, Lₛᵢᴰ, L_Feᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)
    
    μₚ = μ⁰ₘₐₓ*fₚ(T)
    
    Lₗᵢₘ₁ᴰˢⁱ = K_mondo(Si, Kₛᵢ¹)    #23c
    Lₗᵢₘ₂ᴰˢⁱ = ϕ < 0 ? K_mondo((Si)^3, (Kₛᵢ²)^3) : 0    #23d
    
    Fₗᵢₘ₁ᴰˢⁱ = min((μᴰ)/(μₚ*Lₗᵢₘᴰ), Lₚₒ₄ᴰ, Lₙᴰ, L_Feᴰ)  #23a
    Fₗᵢₘ₂ᴰˢⁱ = min(1, 2.2*max(0, Lₗᵢₘ₁ᴰˢⁱ - 0.5)) #23b

    return θₘˢⁱᴰ*Lₗᵢₘ₁ᴰˢⁱ*min(5.4, ((4.4*exp(-4.23*Fₗᵢₘ₁ᴰˢⁱ)*Fₗᵢₘ₂ᴰˢⁱ + 1)*(1 + 2*Lₗᵢₘ₂ᴰˢⁱ)))   #22
end





@inline function (pisces::PISCES)(::Val{:P}, x, y, z, t, P, Z, M, POC, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, L_day, PARᴾ, T, zₘₓₗ, zₑᵤ) 
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`
    δᴾ = bgc.exudation_of_DOC[1]
    mᴾ = bgc.phytoplankton_mortality_rate[1]
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    αᴾ = bgc.initial_slope_of_PI_curve[1]
    
    #equaitons here
    sh = 

    [pₚᶻ, pₚᴹ] = bgc.preference_for_nanophytoplankton
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    gₚᶻ = gᴶ(P, pₚᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    gₚᴹ =  gᴶ(P, pₚᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)

    t_darkᴾ = 

    Lₗᵢₘᴾ, Lₚₒ₄ᴾ, Lₙₕ₄ᴾ, Lₙₒ₃ᴾ, Lₙᴾ, L_Feᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)
    
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ)

    return (1-δᴾ)*μᴾ*P - mᴾ*K_mondo(P, Kₘ)*P - sh*wᴾ*P^2 - gₚᶻ*Z - gₚᴹ*M    #eq 1
end

@inline function (pisces::PISCES)(::Val{:D}, x, y, z, t, D, Z, M, POC, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, L_day, PARᴰ, T, zₘₓₗ, zₑᵤ)
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`
    δᴰ = bgc.exudation_of_DOC[2]
    mᴰ = bgc.phytoplankton_mortality_rate[2]
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms
    αᴰ = bgc.initial_slope_of_PI_curve[2]

    #equaitons here
    sh = 

    [p_Dᶻ, p_Dᴹ] = bgc.preference_for_diatoms
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    g_Dᶻ = gᴶ(D, p_Dᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    g_Dᴹ =  gᴶ(D, p_Dᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)
 
    t_darkᴰ = 

    Lₗᵢₘᴰ, Lₚₒ₄ᴰ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ, Lₙᴰ, Lₛᵢᴰ, L_Feᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)

    wᴰ = wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ) #13
    
    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)

    return (1-δᴰ)*μᴰ*D - mᴰ*K_mondo(D, Kₘ)*D - sh*wᴰ*D^2 - g_Dᶻ*Z - g_Dᴹ*M    #eq 9
end

@inline function (pisces:PISCES)(::Val{:Pᶜʰˡ}, x, y, z, t, P, Z, M, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PARᴾ, T, L_day, zₘₓₗ, zₑᵤ)
    δᴾ = bgc.exudation_of_DOC[1]
    αᴾ = bgc.initial_slope_of_PI_curve[1]
    θₘᵢₙᶜʰˡ = bgc.min_ChlC_ratios_of_phytoplankton
    mᴾ = bgc.phytoplankton_mortality_rate[1]
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton

    sh = 

    [pₚᶻ, pₚᴹ] = bgc.preference_for_nanophytoplankton
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    gₚᶻ = gᴶ(P, pₚᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    gₚᴹ =  gᴶ(P, pₚᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)

    t_darkᴾ = 
    
    Lₗᵢₘᴾ, Lₚₒ₄ᴾ, Lₙₕ₄ᴾ, Lₙₒ₃ᴾ, Lₙᴾ, L_Feᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)

    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ)

    μ̌ᴾ = μᴾ / f₁(L_day) #15b
    ρᴾᶜʰˡ = 144*μ̌ᴾ * P / (αᴾ* Pᶜʰˡ* (PARᴾ)/L_day) #15a

    return (1-δᴾ)*(12*θₘᵢₙᶜʰˡ + (θₘₐₓᶜʰˡᴾ - θₘᵢₙᶜʰˡ)*ρᴾᶜʰˡ)*μᴾ*P - mᴾ*K_mondo(P, Kₘ)*Pᶜʰˡ - sh*wᴾ*P*Pᶜʰˡ - θ(Pᶜʰˡ, P)*gₚᶻ*Z - θ(Pᶜʰˡ, P)*gₚᴹ*M  #14
end

@inline function (pisces:PISCES)(::Val{:Dᶜʰˡ}, x, y, z, t, D, Dᶜʰˡ, Z, M, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, PARᴰ, T, L_day, zₘₓₗ, zₑᵤ)
    δᴾ = bgc.exudation_of_DOC[2]
    αᴾ = bgc.initial_slope_of_PI_curve[2]
    θₘᵢₙᶜʰˡ = bgc.min_ChlC_ratios_of_phytoplankton
    mᴾ = bgc.phytoplankton_mortality_rate[2]
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton

    sh = 

    [p_Dᶻ, p_Dᴹ] = bgc.preference_for_diatoms
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    g_Dᶻ = gᴶ(D, p_Dᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    g_Dᴹ =  gᴶ(D, p_Dᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)
 
    t_darkᴰ = 

    Lₗᵢₘᴰ, Lₚₒ₄ᴰ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ, Lₙᴰ, Lₛᵢᴰ, L_Feᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)

    wᴰ = wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ) #13

    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)

    μ̌ᴰ = μᴰ / f₁(L_day) #15b
    ρᴰᶜʰˡ = 144*μ̌ᴰ * D / (αᴰ* Dᶜʰˡ* (PARᴰ)/L_day) #15a
  
    return (1-δᴰ)*(12*θₘᵢₙᶜʰˡ + (θₘₐₓᶜʰˡᴰ - θₘᵢₙᶜʰˡ)*ρᴰᶜʰˡ)*μᴰ*D - mᴰ*K_mondo(D, Kₘ)*Dᶜʰˡ - sh*wᴰ*D*Dᶜʰˡ - θ(Dᶜʰˡ, D)*g_Dᶻ*Z - θ(Dᶜʰˡ, D)*g_Dᴹ*M    #14
end

@inline function (pisces:PISCES)(::Val{:Pᶠᵉ}, x, y, z, t, P, Z, M, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)
    δᴾ = bgc.exudation_of_DOC[1]
    θₘₐₓᶠᵉᵖ = bgc.max_iron_quota[1]
    mᴾ = bgc.phytoplankton_mortality_rate[1]
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton[1]
    K_Feᴾᶠᵉᵐⁱⁿ = bgc.min_half_saturation_const_for_iron_uptake[1] # this seems wrong as doesn't quite match parameter list

    Lₗᵢₘᴾ, Lₚₒ₄ᴾ, Lₙₕ₄ᴾ, Lₙₒ₃ᴾ, Lₙᴾ, L_Feᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)

    sh = 

    [pₚᶻ, pₚᴹ] = bgc.preference_for_nanophytoplankton
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    gₚᶻ = gᴶ(P, pₚᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    gₚᴹ =  gᴶ(P, pₚᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)

    μᴾᶠᵉ = μᴵᶠᵉ(P, Pᶠᵉ, θₘₐₓᶠᵉᵖ, Sᵣₐₜᴾ, K_Feᴾᶠᵉᵐⁱⁿ, Pₘₐₓ, L_Feᴾ, bFe)

    return (1-δᴾ)*μᴾᶠᵉ*P - mᴾ*K_mondo(P, Kₘ)*Pᶠᵉ - sh*wᴾ*P*Pᶠᵉ - θ(Pᶠᵉ, P)*gₚᶻ*Z - θ(Pᶠᵉ, P)*gₚᴹ*M  #16
end

@inline function (pisces:PISCES)(::Val{:Dᶠᵉ}, x, y, z, t, D, Z, M, PO₄, Si, NO₃, NH₄, Dᶜʰˡ, Dᶠᵉ)
    δᴰ = bgc.exudation_of_DOC[2]
    θₘₐₓᶠᵉᴰ = bgc.max_iron_quota[2]
    mᴰ = bgc.phytoplankton_mortality_rate[2]
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    Sᵣₐₜᴰ = bgc.size_ratio_of_phytoplankton[2]
    Dₘₐₓ = 

    Lₗᵢₘᴰ, Lₚₒ₄ᴰ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ, Lₙᴰ, Lₛᵢᴰ, L_Feᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)

    wᴰ = wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ) #13

    sh = 

    [p_Dᶻ, p_Dᴹ] = bgc.preference_for_diatoms
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    g_Dᶻ = gᴶ(D, p_Dᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    g_Dᴹ =  gᴶ(D, p_Dᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)
   
    μᴰᶠᵉ = μᴵᶠᵉ(D, Dᶠᵉ, θₘₐₓᶠᵉᴰ, Sᵣₐₜᴰ, K_Feᴰᶠᵉᵐⁱⁿ, Dₘₐₓ, L_Feᴰ, bFe)

    return (1-δᴰ)*μᴰᶠᵉ*D - mᴰ*K_mondo(D, Kₘ)*Dᶠᵉ - sh*wᴰ*D*Dᶠᵉ - θ(Dᶠᵉ, D)*g_Dᶻ*Z - θ(Dᶠᵉ, D)*g_Dᴹ*M    #16
end

@inline function (pisces:PISCES)(::Val{:Dˢⁱ}, D, Dˢⁱ, M, Z, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, PARᴰ, L_day, T, ϕ, zₘₓₗ, zₑᵤ)       #ϕ is latitude
    δᴰ = bgc.exudation_of_DOC[2]
    mᴰ = bgc.phytoplankton_mortality_rate[2]
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms
    αᴰ = bgc.initial_slope_of_PI_curve[2]

    sh = 

    [p_Dᶻ, p_Dᴹ] = bgc.preference_for_diatoms
    Jₜₕᵣₑₛₕᶻ = bgc.specific_food_thresholds_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    grazing_arg_m = grazing_argᴹ(P, POC, D, Z, T)
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    g_Dᶻ = gᴶ(D, p_Dᶻ, Jₜₕᵣₑₛₕᶻ, grazing_arg_z)
    g_Dᴹ =  gᴶ(D, p_Dᴹ, Jₜₕᵣₑₛₕᴹ, grazing_arg_m)
 
    t_darkᴰ = 

    θₒₚₜˢⁱᴰ = fθₒₚₜˢⁱᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, μᴰ, T, ϕ)

    Lₗᵢₘᴰ, Lₚₒ₄ᴰ, Lₙₕ₄ᴰ, Lₙₒ₃ᴰ, Lₙᴰ, Lₛᵢᴰ, L_Feᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)

    wᴰ = wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ) #13
    
    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)

    return θₒₚₜˢⁱᴰ*(1-δᴰ)*μᴰ*D - θ(Dˢⁱ, D)*g_Dᴹ*M -  θ(Dˢⁱ, D)*g_Dᶻ*Z - mᴰ*K_mondo(D, Kₘ)*Dˢⁱ - sh*wᴰ*D*Dˢⁱ #21
end
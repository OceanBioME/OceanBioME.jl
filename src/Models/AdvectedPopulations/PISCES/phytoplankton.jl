@inline θ(I,J) = I/(J + eps(0.0))
@inline K_mondo(I, J) = I/(I + J + eps(0.0))
@inline Cₚᵣₒ(I, Iᶜʰˡ, PARᴵ, L_day, αᴵ, μₚ, Lₗᵢₘᴵ)=1-exp(-αᴵ*(θ(Iᶜʰˡ,I))*PARᴵ/(L_day*μₚ*Lₗᵢₘᴵ + eps(0.0)))

@inline f₁(L_day) = 1.5*L_day/(0.5+L_day)  #eq 3a
@inline t_dark(zₘₓₗ, zₑᵤ) = max(0, zₘₓₗ-zₑᵤ)^2/86400#eq 3b,c
@inline f₂(zₘₓₗ, zₑᵤ, t_dark_lim) = 1 - t_dark(zₘₓₗ, zₑᵤ)/(t_dark(zₘₓₗ, zₑᵤ)+t_dark_lim) #eq 3d

@inline fₚ(T) = bgc.temperature_sensitivity_of_growth^T #eq 4a

@inline L_NH₄(NO₃, NH₄, Kₙₒ₃ᴵ, Kₙₕ₄ᴵ) = Kₙₒ₃ᴵ*NH₄/(Kₙₒ₃ᴵ*Kₙₕ₄ᴵ+Kₙₕ₄ᴵ*NO₃+Kₙₒ₃ᴵ*NH₄) #eq 6d
@inline L_NO₃(NO₃, NH₄, Kₙₒ₃ᴵ, Kₙₕ₄ᴵ) = Kₙₕ₄ᴵ*NO₃/(Kₙₒ₃ᴵ*Kₙₕ₄ᴵ+Kₙₕ₄ᴵ*NO₃+Kₙₒ₃ᴵ*NH₄) #eq 6e
@inline L_Fe(P, Pᶠᵉ, θₒₚₜᶠᵉᵖ, θₘᵢₙᶠᵉᵖ) = min(1, max(0, (θ(Pᶠᵉ, P) - θₘᵢₙᶠᵉᵖ)/θₒₚₜᶠᵉᵖ)) #eq 6f

@inline θᶠᵉₘᵢₙ(P, Pᶜʰˡ, Lₙᴾ, Lₙₒ₃ᴾ) = 0.0016/(55.85) * θ(Pᶜʰˡ, P) + 1.21e-5*14*Lₙᴾ/(55.85*7.625)*1.5+1.15*14*Lₙₒ₃ᴾ/(55.85*7.625) #eq 20 -> Lₙ could be meant to be L_NH₄?

@inline I₁(I, Iₘₐₓ) = min(I, Iₘₐₓ) #eq 7a
@inline I₂(I, Iₘₐₓ) = max(0, I - Iₘₐₓ) #eq 7b
@inline Kᵢᴶ(Kᵢᴶᵐⁱⁿ, J₁, J₂, Sᵣₐₜᴶ) = Kᵢᴶᵐⁱⁿ* (J₁ + Sᵣₐₜᴶ* J₂)/(J₁ + J₂) #eq 7c

@inline function μᴵ(I, Iᶜʰˡ, PARᴵ, L_day, T, αᴵ, Lₗᵢₘᴵ)
    
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero


    μₚ = μ⁰ₘₐₓ*fₚ(T)    #eq 4b

    return μₚ * f₁(L_day) * f₂(zₘₓₗ, zₑᵤ, t_dark_lim) * Cₚᵣₒ(I, Iᶜʰˡ, PARᴵ, L_day, αᴵ, μₚ, Lₗᵢₘᴵ) * Lₗᵢₘᴵ #2b 
end

@inline function (pisces::PISCES)(::Val{:P}, x, y, z, t, P, Z, M, PAR) 
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`
    δᴾ = bgc.exudiation_of_DOC[1]
    mᴾ = bgc.phytoplankton_mortality_rate[1]
    Kₘ = bgc.half_saturation_const_for_mortality
    wᴾ = bgc.min_quadratic_mortality_of_phytoplankton
    αᴾ = bgc.initial_slope_of_PI_curve[1]
    θₒₚₜᶠᵉᵖ = bgc.optimal_iron_quota[1]
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton[1]
    Kₙₒ₃ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_nitrate[1]
    Kₙₕ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_ammonium[1]
    #equaitons here
    sh = 
    gₚᶻ = 
    gₚᴹ =  

    P₁ = I₁(P, Pₘₐₓ)
    P₂ = I₂(P, Pₘₐₓ)

    Kₙₒ₃ᴾ = Kᵢᴶ(Kₙₒ₃ᴾᵐⁱⁿ, P₁, P₂, Sᵣₐₜᴾ)
    Kₙₕ₄ᴾ = Kᵢᴶ(Kₙₕ₄ᴾᵐⁱⁿ, P₁, P₂, Sᵣₐₜᴾ)

    Lₚₒ₄ᴾ = K_mondo(PO₄, Kₚₒ₄ᴾ) #6b
    Lₙₕ₄ᴾ = L_NH₄(NO₃, NH₄, Kₙₒ₃ᴾ, Kₙₕ₄ᴾ)
    Lₙₒ₃ᴾ = L_NO₃(NO₃, NH₄, Kₙₒ₃ᴾ, Kₙₕ₄ᴾ)
    Lₙᴾ = Lₙₒ₃ᴾ + Lₙₕ₄ᴾ         #6c

    θₘᵢₙᶠᵉᵖ = θᶠᵉₘᵢₙ(P, Pᶜʰˡ, Lₙᴾ, Lₙₒ₃ᴾ)
    L_Feᴾ = L_Fe(P, Pᶠᵉ ,θₒₚₜᶠᵉᵖ, θₘᵢₙᶠᵉᵖ)

    Lₗᵢₘᴾ = min(Lₚₒ₄ᴾ, Lₙᴾ, L_Feᴾ) #6a
    
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ)

    return (1-δᴾ)*μᴾ*P - mᴾ*K_mondo(P, Kₘ)*P - sh*wᴾ*P^2 - gₚᶻ*Z - gₚᴹ*M    #eq 1
end
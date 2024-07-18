#This document contains functions for:
    #λₚₛᵢ¹ (eq52, parametrisation of dissolution rate of PSi)
    #Forcing for PSi (eq51)

@inline function χ_lab(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z)

    χ_lab⁰ = bgc.proportion_of_the_most_labile_phase_in_PSi
    zₘₐₓ = max(zₘₓₗ, zₑᵤ)

    if z <= zₘₓₗ
        return χ_lab⁰
    else
        return χ_lab⁰*exp(-(λₚₛᵢˡᵃᵇ - λₚₛᵢʳᵉᶠ)*((z-zₘₐₓ)/w_GOC(zₑᵤ, zₘₓₗ))) #eq53
    end
end

@inline function λₚₛᵢ¹(zₘₓₗ, zₑᵤ, z, T, Si)

    λₚₛᵢˡᵃᵇ = bgc.fast_dissolution_rate_of_BSi
    λₚₛᵢʳᵉᶠ = bgc.slow_dissolution_rate_of_BSi

    λₚₛᵢ = χ_lab(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z)*λₚₛᵢˡᵃᵇ + (1 - χ_lab(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z))*λₚₛᵢʳᵉᶠ

    Si_eq = 10^(6.44 - 968/(T + 273.15))
    Siₛₐₜ = (Si_eq - Si)/Si_eq

    return λₚₛᵢ*(0.225*(1 + T/15)*Siₛₐₜ + 0.775*(((1 + T/400)^4)*Siₛₐₜ)^9) #eq52
end

@inline function (pisces::PISCES)(::Val{:PSi}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅, D_dust) 

    Kₘ = bgc.half_saturation_const_for_mortality
    wₘₐₓᴰ = bgc.max_quadratic_mortality_of_diatoms.D
    αᴰ = bgc.initial_slope_of_PI_curve.D
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    Dissₛᵢ =

    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)

    PARᴰ = PARᴰ(PAR¹, PAR², PAR³)

    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅)[1]

    wᴰ = wᴾ + wₘₐₓᴰ*(1-Lₗᵢₘᴰ)

    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)

    θˢⁱᴰ = fθₒₚₜˢⁱᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, μᴰ, T, ϕ, Si̅)

    return  θˢⁱᴰ*grazingᴹ(P, D, Z, POC, T)[3]*M +  θˢⁱᴰ*grazingᶻ(P, D, Z, POC, T)[3]*Z + 
     θˢⁱᴰ*mᴰ*K_mondo(D, Kₘ)*Dˢⁱ + sh*wᴰ*D*Dˢⁱ - λₚₛᵢ¹(zₘₓₗ, zₑᵤ, z, T, Si)*Dissₛᵢ*PSi #add partial derivative here
end
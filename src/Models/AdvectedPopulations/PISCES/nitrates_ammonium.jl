
#This document contains functions for:
    #μₙₒ₃ᴶ,  μₙₕ₄ᴶ (eq8)
    #ΔO₂ (eq57)
    #Nitrif (eq56)
    #N_fix (eq58)
    #Forcing for NO₃ and NH₄ (eqs54, 55)

@inline function μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T,  zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ) 
    αᴾ = bgc.initial_slope_of_PI_curve.P
    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[1]
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ) 
    Lₙₒ₃ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[4]
    Lₙₕ₄ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[3]
    return μᴾ * K_mondo(Lₙₒ₃ᴾ, Lₙₕ₄ᴾ) #eq8
end

@inline function μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ)
    αᴾ = bgc.initial_slope_of_PI_curve.P
    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[1]
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ) 
    Lₙₒ₃ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[4]
    Lₙₕ₄ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[3]
    return μᴾ * K_mondo(Lₙₕ₄ᴾ, Lₙₒ₃ᴾ) #eq8
end

@inline function μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ) 
    αᴰ = bgc.initial_slope_of_PI_curve.D
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[1]
    μᴰ =  μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)
    Lₙₒ₃ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[4]
    Lₙₕ₄ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[3]
    return μᴰ * K_mondo(Lₙₒ₃ᴰ, Lₙₕ₄ᴰ) #eq8
end

@inline function μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ)
    αᴰ = bgc.initial_slope_of_PI_curve[2]
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[1]
    μᴰ =  μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)
    Lₙₒ₃ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[4]
    Lₙₕ₄ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[3]
    return μᴰ * K_mondo(Lₙₕ₄ᴰ, Lₙₒ₃ᴰ) #eq8
end

@inline function ΔO₂(O₂)
    O₂ᵐⁱⁿ¹ = bgc.half_sat_const_for_denitrification1
    O₂ᵐⁱⁿ² = bgc.half_sat_const_for_denitrification2

    return min(1, max(0.4*(O₂ᵐⁱⁿ¹-O₂)/(O₂ᵐⁱⁿ²+O₂))) #eq57
end

@inline Nitrif(NH₄, O₂, λₙₕ₄, PAR) = λₙₕ₄*NH₄*(1-ΔO₂(O₂))/(1+PAR) #eq56a

@inline function (pisces::PISCES)(::Val{:NO₃}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅) 

    λₙₕ₄ =  bgc.max_nitrification_rate

    Rₙₕ₄ = bgc.NC_stoichiometric_ratio_of_ANOTHERPLACEHOLDER
    Rₙₒ₃ = bgc.NC_stoichiometric_ratio_of_dentitrification

    μₙₒ₃ᴾ = μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ)
    μₙₒ₃ᴰ = μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ)

    return Nitrif(NH₄, O₂, λₙₕ₄, PAR) - μₙₒ₃ᴾ*P - μₙₒ₃ᴰ*D - Rₙₕ₄*λₙₕ₄*ΔO₂(O₂)*NH₄ - Rₙₒ₃*Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact)
end

# The following relate specifically to NH₄ forcing
#Change to ifelse

@inline Lₙᴰᶻ(Lₙᴾ) = ifelse(Lₙᴾ>=0.08, 0.01, 1 - Lₙᴾ) #eq58
    
@inline function N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR) #eq 58b
    N_fixᵐ = bgc.max_rate_of_nitrogen_fixation
    K_Feᴰᶻ = bgc.Fe_half_saturation_constant_of_nitrogen_fixation
    Kₚₒ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_phosphate.P
    E_fix = bgc.photosynthetic_parameter_of_nitrogen_fixation
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero
    μₚ = μ⁰ₘₐₓ*fₚ(T)
    Lₙᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[5]

    return N_fixᵐ*max(0,μₚ - 2.15)*Lₙᴰᶻ(Lₙᴾ)*min(K_mondo(bFe, K_Feᴰᶻ), K_mondo(PO₄, Kₚₒ₄ᴾᵐⁱⁿ))*(1 - e^{-PAR/E_fix})
end

@inline function (pisces::PISCES)(::Val{:NH₄}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅) 
    
    #Parameters
    γᶻ = bgc.excretion_as_DOM.Z
    σᶻ = bgc.non_assimilated_fraction.Z
    γᴹ = bgc.excretion_as_DOM.M
    σᴹ = bgc.non_assimilated_fraction.M
    λₙₕ₄ = bgc.max_nitrification_rate
    t_darkᴾ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.P
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D

    bFe = 1 #Change this!

    #L_day
    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)
    
    #Grazing
    grazingᶻ = grazingᶻ(P, D, POC, T)
    grazingᴹ = grazingᴹ(P, D, Z, POC, T)
    ∑gᶻ = grazingᶻ[1]
    ∑gᴹ = grazingᴹ[1]
    ∑g_FFᴹ = ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC)

    #Gross growth efficiency
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z)
    eᴹ =  eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ,Pᶠᵉ, Dᶠᵉ, SFe, P, D, POC)

    μₙₕ₄ᴾ = μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴾ, t_darkᴾ)
    μₙₕ₄ᴰ = μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day, PARᴰ, t_darkᴰ)
   
    return γᶻ*(1-eᶻ-σᶻ)*∑gᶻ*Z + γᴹ*(1-eᴹ-σᴹ)*(∑gᴹ + ∑g_FFᴹ)*M + γᴹ*Rᵤₚᴹ(M, T) + Remin(O₂, NO₃, PO₄, NH₄, DOC, T, bFe, Bact) + Denit(NO₃, PO₄, NH₄, DOC, O₂, T, bFe, Bact) + N_fix(bFe, PO₄, T, P, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, PAR) - Nitrif(NH₄, O₂, λₙₕ₄, PAR) - λₙₕ₄*ΔO₂(O₂)*NH₄ - μₙₕ₄ᴾ*P - μₙₕ₄ᴰ*D
end

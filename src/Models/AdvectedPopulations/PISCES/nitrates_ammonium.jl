 # TO DO:   
    #Fill in parameters for functions
    #Where is Pₘₐₓ defined? 
    #Write PAR̄, where to get PAR̄₁? (56b)
    #Set values for Rₙₕ₄ and Rₙₒ₃.


# For use in NO₃ and NH₄ forcing equations.

@inline function μₙₒ₃ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T,  zₘₓₗ, zₑᵤ, L_day) 
    #PARᴾ = 
    t_darkᴾ = 
    αᴾ = bgc.initial_slope_of_PI_curve[1]
    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[1]
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ) 
    Lₙₒ₃ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[4]
    Lₙₕ₄ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[3]
    return μᴾ * K_mondo(Lₙₒ₃ᴾ, Lₙₕ₄ᴾ) #eq8
end

@inline function μₙₕ₄ᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day)
    #PARᴾ = 
    t_darkᴾ = 
    αᴾ = bgc.initial_slope_of_PI_curve[1]
    Lₗᵢₘᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[1]
    μᴾ = μᴵ(P, Pᶜʰˡ, PARᴾ, L_day, T, αᴾ, Lₗᵢₘᴾ, zₘₓₗ, zₑᵤ, t_darkᴾ) 
    Lₙₒ₃ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[4]
    Lₙₕ₄ᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[3]
    return μᴾ * K_mondo(Lₙₕ₄ᴾ, Lₙₒ₃ᴾ) #eq8
end

@inline function μₙₒ₃ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day) 
    #PARᴰ = 
    t_darkᴰ = 
    αᴰ = bgc.initial_slope_of_PI_curve[2]
    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[1]
    μᴰ =  μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)
    Lₙₒ₃ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[4]
    Lₙₕ₄ᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[3]
    return μᴰ * K_mondo(Lₙₒ₃ᴰ, Lₙₕ₄ᴰ) #eq8
end

@inline function μₙₕ₄ᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, T, zₘₓₗ, zₑᵤ, L_day)
    #PARᴰ = 
    t_darkᴰ = 
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

PAR̄() = 0 #eq56b

@inline function Nitrif(NH₄, O₂, λₙₕ₄) 
    O₂ᵐⁱⁿ¹ = bgc.half_sat_const_for_denitrification1
    O₂ᵐⁱⁿ² = bgc.half_sat_const_for_denitrification2
    
    return λₙₕ₄*NH₄*(1-ΔO₂(O₂))/(1+PAR̄()) #eq56a

# For NO₃ forcing only
Rₙₕ₄ = 0 # set this value
Rₙₒ₃ = 0.86 # check this value

@inline function (pisces::PISCES)(::Val{:NO₃}, x, y, z, t, P, D, NH₄, O₂, PAR) 

    λₙₕ₄ =  bgc.max_nitrification_rate

    return Nitrif(λₙₕ₄, NH₄, O₂) - μₙₒ₃ᴾ()*P - μₙₒ₃ᴰ()*D - Rₙₕ₄*λₙₕ₄*ΔO₂(O₂)*NH₄ - Rₙₒ₃*Denit()
end

# The following relate specifically to NH₄ forcing

@inline function Lₙᴰᶻ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ) #eq58a
    Lₙᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ)[5]
    #Lₗᵢₘᴾ, Lₚₒ₄ᴾ, Lₙₕ₄ᴾ, Lₙₒ₃ᴾ, Lₙᴾ, L_Feᴾ = Lᴾ(P, PO₄, NO₃, NH₄, Pᶜʰˡ, Pᶠᵉ), check the correct way to call the function.
    if Lₙᴾ >= 0.08 #Fill parameters for this function. Check where Pₘₐₓ defined.
        return 0.01
    else
        return 1 - Lₙᴾ
    end
        
@inline function N_fix(bFe, PO₄, PAR, T) #eq 58b
    N_fixᵐ = bgc.max_rate_of_nitrogen_fixation
    K_Feᴰᶻ = bgc.Fe_half_saturation_constant_of_nitrogen_fixation
    Kₚₒ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_phosphate[1]
    E_fix = bgc.photosynthetic_parameter_of_nitrogen_fixation
    μ⁰ₘₐₓ = bgc.growth_rate_at_zero
    μₚ = μ⁰ₘₐₓ*fₚ(T)

    return N_fixᵐ*max(0,μₚ - 2.15)*Lₙᴰᶻ()*min(K_mondo(bFe, K_Feᴰᶻ), K_mondo(PO₄, Kₚₒ₄ᴾᵐⁱⁿ))*(1 - e^{-PAR/E_fix})
end


@inline function (pisces::PISCES)(::Val{:NH₄}, x, y, z, t, P, D, NH₄, O₂, bFe, POC, GOC, PAR) 
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`

    γᶻ = bgc.excretion_as_DOM[1]
    σᶻ = bgc.non_assimilated_fraction[1]
    γᴹ = bgc.excretion_as_DOM[2]
    σᴹ = bgc.non_assimilated_fraction[2]
    λₙₕ₄ = bgc.max_nitrification_rate

    pₚᶻ = bgc.preference_for_nanophytoplankton[1]
    p_Dᶻ = bgc.preference_for_diatoms[1]
    p_pocᶻ = bgc.preference_for_POC[1]
    pₚᴹ = bgc.preference_for_nanophytoplankton[2]
    p_Dᴹ = bgc.preference_for_diatoms[2]
    pₚₒᴹ = bgc.preference_for_POC[2]
    p_zᴹ = bgc.preference_for_microzooplankton
    Jₜₕᵣₑₛₕᴹ = bgc.specific_food_thresholds_for_mesozooplankton
    grazing_arg_z = grazing_argᶻ(P, POC, D, T) 
    grazing_arg_m = grazing_argᴹ(P, POC, D, T) 
    
    grazingᶻ = grazingᶻ()

    gₚᶻ = grazingᶻ[2]
    g_Dᶻ = grazingᶻ[3]
    gₚₒᶻ = grazingᶻ[4]
    g_Zᴹ = grazingᴹ[5]
    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z)

    gₚᴹ = grazingᴹ[2]
    g_Dᴹ = grazingᴹ)[3]
    gₚₒᴹ = grazingᴹ[4]
    eᴹ = eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ, N, Fe, P, D, POC, Z)

    ∑gᶻ = grazingᶻ[1]
    ∑gᴹ = grazingᴹ[1]
   
    return γᶻ*(1-eᶻ-σᶻ)*∑gᶻ*Z + γᴹ*(1-eᴹ-σᴹ)*(∑gᴹ + ∑g_FFᴹ())*M + γᴹ*Rᵤₚᴹ() + Remin() + Denit() + N_fix() - Nitrif() - λₙₕ₄*ΔO₂()*NH₄ - μₙₕ₄ᴾ()*P - μₙₕ₄ᴰ()*D
end

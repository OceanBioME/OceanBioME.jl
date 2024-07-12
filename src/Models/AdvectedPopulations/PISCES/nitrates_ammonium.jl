 # TO DO:   
    #Fill in the variables for these functions once they are written. Add relevant parameters to parameter lists inside forcing function. 
    #Fill arguments for μᴾₙₒ₃(), μᴰₙₒ₃() (8), Denit() (33b) - for NO₃.
    #Fill arguments for Lₙᴾ() (6c), μₚ() (4b), gᴹ() (26a), eᶻ() (27), Lₙᴰᶻ() (58a) - for NH₄.
    #Where is Pₘₐₓ defined? 
    #Write PAR̄, where to get PAR̄₁? (56b)
    #Set values for Rₙₕ₄ and Rₙₒ₃.


# For use in NO₃ and NH₄ forcing equations.

@inline ΔO₂(O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂) = min(1, max(0.4*(O₂ᵐⁱⁿ¹-O₂)/(O₂ᵐⁱⁿ²+O₂))) #(57)
PAR̄() = 0 #(56b)

@inline Nitrif(λₙₕ₄, NH₄, O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂) = λₙₕ₄*NH₄*(1-ΔO₂(O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂))/(1+PAR̄()) #(56a)

# For NO₃ forcing only
Rₙₕ₄ = 0 # set this value
Rₙₒ₃ = 0.86 # check this value

#NO₃ forcing also requires μᴾₙₒ₃, μᴰₙₒ₃ (8), Denit (33b). 

@inline function (pisces::PISCES)(::Val{:NO₃}, x, y, z, t, P, D, NH₄, O₂, PAR) 

    O₂ᵐⁱⁿ¹ = bgc.half_sat_const_for_denitrification1
    O₂ᵐⁱⁿ² = bgc.half_sat_const_for_denitrification2
    λₙₕ₄ =  bgc.max_nitrification_rate

    return Nitrif(λₙₕ₄, NH₄, O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂) - μₙₒ₃ᴾ()*P - μₙₒ₃ᴰ()*D - Rₙₕ₄*λₙₕ₄*ΔO₂(O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂)*NH₄ - Rₙₒ₃*Denit()
end


# The following relate specifically to NH₄ forcing

@kwdef function Lₙᴰᶻ() #(58a), Lₙᵖ (6).
    if Lₙᴾ() >= 0.08 #Fill parameters for this function. Check where Pₘₐₓ defined.
        return 0.01
    else
        return 1 - Lₙᴾ()
        
@inline N_fix(N_fixᵐ, K_Feᴰᶻ, Kₚₒ₄ᴾᵐⁱⁿ, E_fix, Pₘₐₓ, Kᵢᴾᵐⁱⁿ, Sᵣₐₜᴾ) = N_fixᵐ*max(0,μₚ() - 2.15)*Lₙᴰᶻ()*min(bFe/(K_Feᴰᶻ + bFe), PO₄/(Kₚₒ₄ᴾᵐⁱⁿ + PO₄))*(1 - e^{-PAR/E_fix}) #(58b)


# Define sum of grazing rates, as this quantity freqeuently appears

@inline ∑gᴹ() = gᴹ(P, ) + gᴹ(D, ) + gᴹ(Z, ) + gᴹ(POC, )

# NH₄ forcing also requires eᶻ, eᴹ (27), Rᵤₚᴹ (30b), gᶻ (26a), g_FF (29), μᴾₙₕ₄, μᴰₙₕ₄() (8)

@inline function (pisces::PISCES)(::Val{:NH₄}, x, y, z, t, P, D, NH₄, O₂, bFe, POC, GOC, PAR) 
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`

    γᶻ = bgc.excretion_as_DOM[1]
    σᶻ = bgc.non_assimilated_fraction[1]
    γᴹ = bgc.excretion_as_DOM[2]
    σᴹ = bgc.non_assimilated_fraction[2]
    λₙₕ₄ = bgc.max_nitrification_rate
    #Required for (57) and (56)
    O₂ᵐⁱⁿ¹ = bgc.half_sat_const_for_denitrification1
    O₂ᵐⁱⁿ² = bgc.half_sat_const_for_denitrification2
    #Required for (58)
    N_fixᵐ = bgc.max_rate_of_nitrogen_fixation
    K_Feᴰᶻ = bgc.Fe_half_saturation_constant_of_nitrogen_fixation
    Kₚₒ₄ᴾᵐⁱⁿ = bgc.min_half_saturation_const_for_phosphate[1]
    E_fix = bgc.photosynthetic_parameter_of_nitrogen_fixation
    Sᵣₐₜᴾ = bgc.size_ratio_of_phytoplankton[1]
     # Pₘₐₓ = ?
  
    return γᶻ*(1-eᶻ()-σᶻ)*∑gᶻ()*Z + γᴹ*(1-eᴹ()-σᴹ)*(∑gᴹ() + g_FFᴹ(POC, ) + g_FFᴹ(GOC, ))*M + γᴹ*Rᵤₚᴹ() + Remin() + Denit() + N_fix(N_fixᵐ, K_Feᴰ, Kₚₒ₄ᴾᵐⁱⁿ, E_fix, Pₘₐₓ, Kᵢᴾᵐⁱⁿ, Sᵣₐₜᴾ) - Nitrif(λₙₕ₄, NH₄, O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂) - λₙₕ₄*ΔO₂(O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂)*NH₄ - μₙₕ₄ᴾ()*P - μₙₕ₄ᴰ()*D
# TO DO: Forcing equations depend on functions from other documents. Fill in the variables for these functions once they are written. Add relevant parameters to parameter lists inside forcing function. Where is Pₘₐₓ defined?


# For use in NO₃ and NH₄ forcing equations.

@inline ΔO₂(O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂) = min(1, max(0.4*(O₂ᵐⁱⁿ¹-O₂)/(O₂ᵐⁱⁿ²+O₂))) #(57)
PAR̄ = 0 # where are PAR₁, etc defined? (56b)

@inline Nitrif(λₙₕ₄, NH₄, O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂) = λₙₕ₄*NH₄*(1-ΔO₂)/(1+PAR̄) #(56a)

# For NO₃ forcing only
Rₙₕ₄ = 0 # set this value
Rₙₒ₃ = 0.86 # check this value

#NO₃ forcing also requires μᴾₙₒ₃, μᴰₙₒ₃ (8), Denit (33b)

@inline function (pisces::PISCES)(::Val{:NO₃}, x, y, z, t, P, D, NH₄, O₂, PAR) #May need to add other tracers here to call with functions defined elsewhere
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`

    #Write parameters as symbols, also check which parameters are required to pass into functions.
    O₂ᵐⁱⁿ¹ = 
    O₂ᵐⁱⁿ² = 

    return Nitrif(λₙₕ₄, NH₄, O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂) - μᴾₙₒ₃()*P - μᴰₙₒ₃()*D - Rₙₕ₄*λₙₕ₄*ΔO₂(O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂)*NH₄ - Rₙₒ₃*Denit()
end


# The following relate specifically to NH₄ forcing

@kwdef function Lₙᴰᶻ(Pₘₐₓ, Kᵢᴾᵐⁱⁿ, Sᵣₐₜᴾ) #(58a), Lₙᵖ (6).
    if Lₙᴾ() >= 0.08 #Fill parameters for this function. Check where Pₘₐₓ defined.
        return 0.01
    else
        return 1 - Lₙᴾ()
        
@inline N_fix(N_fixᵐ, K_Feᴰ, Kₚₒ₄ᴾᵐⁱⁿ, E_fix, Pₘₐₓ, Kᵢᴾᵐⁱⁿ, Sᵣₐₜᴾ) = N_fixᵐ*max(0,μₚ() - 2.15)*Lₙᴰᶻ(Pₘₐₓ, Kᵢᴾᵐⁱⁿ, Sᵣₐₜᴾ)*min(bFe/(K_Feᴰᶻ + bFe), PO₄/(Kₚₒ₄ᴾᵐⁱⁿ + PO₄))*(1 - e^{-PAR/E_fix}) #(58b)


# Define sum of grazing rates, as this quantity freqeuently appears

@inline ∑gᴹ() 

# NH₄ forcing also requires eᶻ, eᴹ (27), Rᵤₚᴹ (30b), gᶻ (26a), g_FF (29), μᴾₙₕ₄, μᴰₙₕ₄() (8)

@inline function (pisces::PISCES)(::Val{:NH₄}, x, y, z, t, P, D, NH₄, O₂, bFe, POC, GOC, PAR) 
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`

    γᶻ = 
    σᶻ = 
    γᴹ = 
    σᴹ = 
    λₙₕ₄ = 
    #Required for (57) and (56)
    O₂ᵐⁱⁿ¹ = 
    O₂ᵐⁱⁿ² =
    #Required for (58)
    N_fixᵐ = 
    K_Feᴰ = 
    Kₚₒ₄ᴾᵐⁱⁿ = 
    E_fix, Pₘₐₓ = 
    Kᵢᴾᵐⁱⁿ = 
    Sᵣₐₜᴾ = 
  
    return γᶻ*(1-eᶻ()-σᶻ)*∑gᴹ()*Z + γᴹ*(1-eᴹ()-σᴹ)*(∑gᴹ() + g_FFᴹ() + g_FFᴹ())*M + γᴹ*Rᵤₚᴹ() + Remin() + Denit() + N_fix(N_fixᵐ, K_Feᴰ, Kₚₒ₄ᴾᵐⁱⁿ, E_fix, Pₘₐₓ, Kᵢᴾᵐⁱⁿ, Sᵣₐₜᴾ) - Nitrif(λₙₕ₄, NH₄, O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂) - λₙₕ₄*ΔO₂(O₂ᵐⁱⁿ¹, O₂ᵐⁱⁿ², O₂)*NH₄ - μₙₕ₄ᴾ()*P - μₙₕ₄ᴰ()*D
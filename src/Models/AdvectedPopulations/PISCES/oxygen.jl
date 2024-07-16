#TO DO:
    #
    
#This document contains functions for:
    #O₂ forcing (eq83)

    @inline function (pisces::PISCES)(::Val{:O₂}, x, y, z, t, P, D, Z, M, PAR) 

        O₂ᵘᵗ = bgc.OC_for_ammonium_based_processes
        O₂ⁿⁱᵗ = bgc.OC_ratio_of_nitrification
        γᶻ = bgc.excretion_as_DOM[1]
        γᴹ = bgc.excretion_as_DOM[2]
        σᶻ = bgc.non_assimilated_fraction[1]
        σᴹ = bgc.non_assimilated_fraction[2]
    
        return O₂ᵘᵗ*(μₙₕ₄ᴾ()*P + μₙₕ₄ᴰ()*D) + (O₂ᵘᵗ + O₂ⁿⁱᵗ)*(μₙₒ₃ᴾ()*P + μₙₒ₃ᴰ()*D) + O₂ⁿⁱᵗ*N_fix() - O₂ᵘᵗ*γᶻ*(1 - eᶻ() - σᶻ)*(∑gᶻ())*Z - O₂ᵘᵗ*γᴹ*(1 - eᴹ() - σᴹ)*(∑gᴹ() + g_FF(POC, ) + g_FF(GOC, ))*M - O₂ᵘᵗ*γᴹ*Rᵤₚᴹ() - O₂ᵘᵗ*Remin() - O₂ⁿⁱᵗ*Nitrif()
    end
#This document contains functions for:
    #Forcing for DIC.
    #Forcing for Alk.

@inline function (pisces::PISCES)(::Val{:DIC}, x, y, z, t, P, D, Z, M, T, POC, GOC, PAR)
    zₑᵤ = 
    zₘₓₗ = 

    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z)
    eᴹ = eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ, N, Fe, P, D, POC, Z)

    ∑gᶻ = grazingᶻ()[1]
    ∑gᴹ = grazingᴹ()[1]

    return γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z + γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC))*M + γᴹ*Rᵤₚᴹ + Remin() + Denit() + λ_CaCO₃¹()*CaCO₃ - P_CaCO₃() - μᴰ()*D - μᴾ()*P #eq59
end


@inline function (pisces::PISCES)(::Val{:Alk}, x, y, z, t, P, D, Z, M, POC, GOC, PAR) # eq59
    
    θᴺᶜ = bgc.NC_redfield_ratio
    rₙₒ₃¹ = bgc. CN_ratio_of_denitrification
    rₙₕ₄¹ = bgc.CN_ratio_of_ammonification
    γᶻ = bgc.excretion_as_DOM[1]
    σᶻ = bgc.non_assimilated_fraction[1]
    γᴹ = bgc.excretion_as_DOM[2]
    σᴹ = bgc.non_assimilated_fraction[2]
    λₙₕ₄ = bgc.max_nitrification_rate

    eᶻ = eᴶ(eₘₐₓᶻ, σᶻ, gₚᶻ, g_Dᶻ, gₚₒᶻ, g_zᴹ, N, Fe, P, D, POC, Z)
    eᴹ = eᴶ(eₘₐₓᴹ, σᴹ, gₚᴹ, g_Dᴹ, gₚₒᴹ, g_zᴹ, N, Fe, P, D, POC, Z)


    ∑gᶻ = grazingᶻ()[1]
    ∑gᴹ = grazingᴹ()[1]

    μₙₒ₃ᴾ = 
    μₙₕ₄ᴾ = 
    μₙₒ₃ᴰ = 
    μₙₕ₄ᴰ = 

   
    return θᴺᶜ*Remin() + θᴺᶜ*(rₙₒ₃¹ + 1)*Denit() + θᴺᶜ*γᶻ*(1 - eᶻ - σᶻ)*∑gᶻ*Z 
    + θᴺᶜ*γᴹ*(1 - eᴹ - σᴹ)*(∑gᴹ + ∑g_FFᴹ(zₑᵤ, zₘₓₗ, T, POC, GOC) + θᴺᶜ*γᴹ*Rᵤₚ())*M
    + θᴺᶜ*μₙₒ₃ᴾ*P + θᴺᶜ*μₙₒ₃ᴰ*D + θᴺᶜ*N_fix() + 2*λ_CaCO₃¹()*CaCO₃
    + θᴺᶜ*ΔO₂()*(rₙₕ₄¹ - 1)*λₙₕ₄*NH₄ - θᴺᶜ*μₙₕ₄ᴾ*P - θᴺᶜ*μₙₕ₄ᴰ*D- 2*θᴺᶜ*Nitrif() - 2*P_CaCO₃
end
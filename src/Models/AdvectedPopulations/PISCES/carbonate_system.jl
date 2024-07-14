#TO DO:
    #Do we assume θᴺᶜ = N/C?
    #Check all of this for typos.


@inline function (pisces::PISCES)(::Val{:DIC}, x, y, z, t, P, D, Z, M, POC, GOC, PAR) # eq59
    


    return γᶻ*(1 - eᶻ() - σᶻ)*∑gᶻ()*Z + γᴹ*(1 - eᴹ() - σᴹ)*(∑gᴹ() + g_FFᴹ(POC, ) + g_FFᴹ(GOC, ))*M
    + γᴹ*Rᵤₚᴹ + Remin() + Denit() + λ_CaCO₃¹*CaCO₃ - P_CaCO₃ - μᴰ()*D - μᴾ()*P
end




@inline function (pisces::PISCES)(::Val{:Alk}, x, y, z, t, P, D, Z, M, POC, GOC, PAR) # eq59
    
    θᴺᶜ = 

    return θᴺᶜ*Remin() + θᴺᶜ*(rₙₒ₃¹ + 1)*Denit() + θᴺᶜ*γᶻ*(1 - eᶻ() - σᶻ)*∑gᶻ()*Z 
    + θᴺᶜ*γᴹ*(1 - eᴹ() - σᶻ)*(∑gᴹ() + g_FFᴹ(POC, ) + g_FF(GOC, ) + θᴺᶜ*γᴹ*Rᵤₚ())*M
    + θᴺᶜ*μₙₒ₃ᴾ*P + θᴺᶜ*μₙₒ₃ᴰ*D + θᴺᶜ*N_fix() + 2*λ_CaCO₃¹*CaCO₃
    + θᴺᶜ*ΔO₂()*(rₙₕ₄¹ - 1)*λₙₕ₄*NH₄ - θᴺᶜ*μₙₕ₄ᴾ*P - θᴺᶜ*μₙₕ₄ᴰ*D- 2*θᴺᶜ*Nitrif() - 2*P_CaCO₃
end
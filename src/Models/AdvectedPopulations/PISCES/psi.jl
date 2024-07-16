
#TO DO:
    #Add partial derivative to eq51
    #What is Dissₛᵢ?
    #What are  λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ? λₚₛᵢˢˡᵒʷ, λₚₛᵢᶠᵃˢᵗ given in parameter list but not used?

#This document contains functions for:
    #λₚₛᵢ¹ (eq52, parametrisation of dissolution rate of PSi)
    #Forcing for PSi (eq51)

@inline function χ_lab(zₘₓₗ, zₑᵤ, χ_lab⁰, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z)
    zₘₐₓ = max(zₘₓₗ, zₑᵤ)
    if z <= zₘₓₗ
        return χ_lab⁰
    else
        return χ_lab⁰*exp(-(λₚₛᵢˡᵃᵇ - λₚₛᵢʳᵉᶠ)*((z-zₘₐₓ)/ω_GOC(zₑᵤ, zₘₓₗ))) #eq53
end

@inline function λₚₛᵢ¹(zₘₓₗ, zₑᵤ, χ_lab⁰, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z, T, Si)
    
    
    λₚₛᵢ = χ_lab()*λₚₛᵢˡᵃᵇ + (1 - χ_lab())*λₚₛᵢʳᵉᶠ

    Si_eq = 10^(6.44 - 968/(T + 273.15))
    Siₛₐₜ = (Si_eq - Si)/Si_eq

    return λₚₛᵢ*(0.225*(1 + T/15)*Siₛₐₜ + 0.775*(((1 + T/400)^4)*Siₛₐₜ)^9) #eq52
end

@inline function (pisces::PISCES)(::Val{:PSi}, x, y, z, t, P, D, Z, M, Dˢⁱ, PSi, PAR) 

    Kₘ = bgc.half_saturation_const_for_mortality
    ωᴰ = #fill this 

    zₑᵤ =  
    zₘₓₗ = 

    θˢⁱᴰ = fθₒₚₜˢⁱᴰ()

    return  θˢⁱᴰ*grazingᴹ()[3]*M +  θˢⁱᴰ*grazingᶻ()[3]*Z +  θˢⁱᴰ*mᴰ*K_mondo(D, Kₘ)*Dˢⁱ + sh*ωᴰ*D*Dˢⁱ - λₚₛᵢ¹()*Dissₛᵢ*PSi -  ω_GOC(zₑᵤ, zₘₓₗ)* #add partial derivative here
end
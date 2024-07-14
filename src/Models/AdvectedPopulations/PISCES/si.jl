#To Do:
    #What is λₚₛᵢ¹, λₚₛᵢ* in notes?
    #What is Dissₛᵢ?
    #Is θₒₚₜᴰˢⁱ

@inline function (pisces::PISCES)(::Val{:Si}, x, y, z, t, D, PSi, PAR) 
   
    δᴰ = bgc.exudation_of_DOC
    
    return λₚₛᵢ¹*Dissₛᵢ*PSi - θₒₚₜˢⁱᴰ()*(1-δᴰ)*μᴰ()*D 
#To Do:
    #What is Dissₛᵢ?

@inline function (pisces::PISCES)(::Val{:Si}, x, y, z, t, D, PSi, Dᶜʰˡ, T, PAR) 
   
    δᴰ = bgc.exudation_of_DOC
    αᴰ = bgc.initial_slope_of_PI_curve[2]

    zₘₓₗ = 
    zₑᵤ = 

    PARᴰ = 
    t_darkᴰ = 
    L_day = 

    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)

    λₚₛᵢ¹ = λₚₛᵢ¹()
    
    return λₚₛᵢ¹*Dissₛᵢ*PSi - fθₒₚₜˢⁱᴰ()*(1-δᴰ)*μᴰ*D 
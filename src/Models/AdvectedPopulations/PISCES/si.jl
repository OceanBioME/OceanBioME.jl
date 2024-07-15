#To Do:
    #What is λₚₛᵢ¹, λₚₛᵢ* in notes?
    #What is Dissₛᵢ?

@inline function (pisces::PISCES)(::Val{:Si}, x, y, z, t, D, PSi, Dᶜʰˡ, T, PAR) 
   
    δᴰ = bgc.exudation_of_DOC
    αᴰ = bgc.initial_slope_of_PI_curve[2]

    zₘₓₗ = 
    zₑᵤ = 

    PARᴰ = 
    t_darkᴰ = 
    L_day = 

    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ)[1]

    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)
    
    return λₚₛᵢ¹*Dissₛᵢ*PSi - fθₒₚₜˢⁱᴰ()*(1-δᴰ)*μᴰ*D 
#To Do:
    #What is Dissₛᵢ?

@inline function (pisces::PISCES)(::Val{:Si}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, O₂, T, PAR, PAR¹, PAR², PAR³, zₘₓₗ, zₑᵤ, Si̅) 
   
    δᴰ = bgc.exudation_of_DOC.D
    αᴰ = bgc.initial_slope_of_PI_curve.D
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    Dissₛᵢ = bgc.

    PARᴰ = PARᴰ(PAR¹, PAR², PAR³)
    
    ϕ₀ = bgc.latitude
    L_day_param = bgc.length_of_day
    ϕ = get_ϕ(ϕ₀, y)
    L_day = get_L_day(ϕ, t, L_day_param)

    Lₗᵢₘᴰ = Lᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅)[1]

    μᴰ = μᴵ(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ)

    λₚₛᵢ¹ = λₚₛᵢ¹(zₘₓₗ, zₑᵤ, z, T, Si)
    
    return λₚₛᵢ¹*Dissₛᵢ*PSi - fθₒₚₜˢⁱᴰ(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, μᴰ, T, ϕ, Si̅)*(1-δᴰ)*μᴰ*D 
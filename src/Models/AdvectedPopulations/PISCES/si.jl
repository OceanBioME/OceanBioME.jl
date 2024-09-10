#The silicon compartment of the model is composed of Dˢⁱ, Si, PSi. Silicon is a limiting nutrient for diatoms, but not phytoplankton.
#Silicon is conserved in the model.

# This documentation contains functions for:
    #Si (eq74)

@inline function (bgc::PISCES)(::Val{:Si}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, PAR, PAR₁, PAR₂, PAR₃) #eq74
    #Parameters
    δᴰ = bgc.exudation_of_DOC.D
    αᴰ = bgc.initial_slope_of_PI_curve.D
    t_darkᴰ = bgc.mean_residence_time_of_phytoplankton_in_unlit_mixed_layer.D
    Dissₛᵢ = bgc.dissolution_rate_of_silicon

    λₚₛᵢ¹ = PSi_dissolution_rate(zₘₓₗ, zₑᵤ, z, T, Si, bgc)
    
    #L_day
    φ = bgc.latitude
    φ = latitude(φ, y)


    L_day = day_length(ϕ, t)
    #Diatom growth
    Lₗᵢₘᴰ = D_nutrient_limitation(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)[1]
    PARᴰ = D_PAR(PAR₁, PAR₂, PAR₃, bgc)
    μᴰ = phytoplankton_growth_rate(D, Dᶜʰˡ, PARᴰ, L_day, T, αᴰ, Lₗᵢₘᴰ, zₘₓₗ, zₑᵤ, t_darkᴰ, bgc)
    
    return λₚₛᵢ¹*Dissₛᵢ*PSi - variation_in_SiC_ratio(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, μᴰ, T, ϕ, Si̅, bgc)*(1-δᴰ)*μᴰ*D #eq74
end
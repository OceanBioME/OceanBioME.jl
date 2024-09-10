
#This document contains functions for:
    #λₚₛᵢ¹ (eq52, parametrisation of dissolution rate of PSi) 
    #Forcing for PSi (eq51)

#Diatom frustule made from 2 silica phases, which dissolve at different rates. We formulate the proportion of fastest dissolving pahse as a function of depth.
@inline function labile_phase(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z, bgc)
    χ_lab⁰ = bgc.proportion_of_the_most_labile_phase_in_PSi
    zₘₐₓ = max(abs(zₘₓₗ), abs(zₑᵤ))

    if abs(z) <= zₘₐₓ
        return χ_lab⁰
    else
        return χ_lab⁰*exp(-(λₚₛᵢˡᵃᵇ - λₚₛᵢʳᵉᶠ)*((abs(z)-zₘₐₓ)/(sinking_speed_of_GOC(z, zₑᵤ, zₘₓₗ, bgc) + eps(0.0)))) #eq53
    end
end

#PSi dissolves to Si. Dissolution rate has following formulation. This is a function of Si and T.
@inline function PSi_dissolution_rate(zₘₓₗ, zₑᵤ, z, T, Si, bgc)
    λₚₛᵢˡᵃᵇ = bgc.fast_dissolution_rate_of_PSi #Discrepancy in labelling. Assumed these are  λₚₛᵢᶠᵃˢᵗ,  λₚₛᵢˢˡᵒʷ from parameter list.
    λₚₛᵢʳᵉᶠ = bgc.slow_dissolution_rate_of_PSi
    
    Si_eq = 10^(6.44 - 968/(T + 273.15)) #eq52
    Siₛₐₜ = (Si_eq - Si)/(Si_eq + eps(0.0))
    λₚₛᵢ = labile_phase(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z, bgc)*λₚₛᵢˡᵃᵇ + (1 - labile_phase(zₘₓₗ, zₑᵤ, λₚₛᵢˡᵃᵇ, λₚₛᵢʳᵉᶠ, z, bgc))*λₚₛᵢʳᵉᶠ #eq53b

    return λₚₛᵢ*(0.225*(1 + T/15)*Siₛₐₜ + 0.775*(((1 + T/400)^4)*Siₛₐₜ)^9) 
end

#Forcing for PSi
@inline function (bgc::PISCES)(::Val{:PSi}, x, y, z, t, P, D, Z, M, Pᶜʰˡ, Dᶜʰˡ, Pᶠᵉ, Dᶠᵉ, Dˢⁱ, DOC, POC, GOC, SFe, BFe, PSi, NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, O₂, T, zₘₓₗ, zₑᵤ, Si̅, D_dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃) 
    #Parameters
    Kₘ = bgc.half_saturation_const_for_mortality
    Dissₛᵢ = bgc.dissolution_rate_of_silicon
    mᴰ = bgc.phytoplankton_mortality_rate.D

    #Also required
    τ₀ = bgc.background_shear
    τₘₓₗ = bgc.mixed_layer_shear

    sh = shear(z, zₘₓₗ, τ₀, τₘₓₗ)
    wᴰ = D_quadratic_mortality(D, PO₄, NO₃, NH₄, Si, Dᶜʰˡ, Dᶠᵉ, Si̅, bgc)
    θˢⁱᴰ = nutrient_quota(Dˢⁱ, D)

    return  (θˢⁱᴰ*grazing_M(P, D, Z, POC, T, bgc)[3]*M +  θˢⁱᴰ*grazing_Z(P, D, POC, T, bgc)[3]*Z 
           + mᴰ*concentration_limitation(D, Kₘ)*Dˢⁱ + sh*wᴰ*D*Dˢⁱ - PSi_dissolution_rate(zₘₓₗ, zₑᵤ, z, T, Si, bgc)*Dissₛᵢ*PSi) #removed θˢⁱᴰ from third term, to conserve silicon
end
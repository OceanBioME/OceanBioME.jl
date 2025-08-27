@inline function (bgc::LOBSTER)(::Val{:sPOC}, x, y, z, t, NO₃, NH₄, Fe, P, Z, sPON, bPON, DON, sPOC, bPOC, DOC, PAR)
    aᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    fᵈ = bgc.fast_sinking_mortality_fraction # really dumb definitions
    fᶻ = bgc.slow_sinking_mortality_fraction
    mᵖ = bgc.phytoplankton_mortality
    μᵈ = bgc.small_detritus_remineralisation_rate
    Rᵖ = bgc.phytoplankton_redfield

    return (((1 - fᵈ) * (1 - aᶻ) * (Gᵖ(P, Z, sPON, gᶻ, p̃, kᶻ) + Gᵈ(P, Z, sPON, gᶻ, p̃, kᶻ))
            + (1 - fᵈ) * mᵖ * P^2
            - Gᵈ(P, Z, sPON, gᶻ, p̃, kᶻ)
            + fᶻ * mᶻ * Z^2) * Rᵖ
            - μᵈ * sPOC)
end

@inline function (bgc::LOBSTER)(::Val{:bPOC}, x, y, z, t, NO₃, NH₄, Fe, P, Z, sPON, bPON, DON, sPOC, bPOC, DOC, PAR)
    aᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    fᵈ = bgc.fast_sinking_mortality_fraction # really dumb definitions
    fᶻ = bgc.slow_sinking_mortality_fraction
    mᵖ = bgc.phytoplankton_mortality
    μᵈᵈ = bgc.large_detritus_remineralisation_rate
    Rᵖ = bgc.phytoplankton_redfield
    η = bgc.zooplankton_calcite_dissolution
    ρᶜᵃᶜᵒ³ = bgc.organic_carbon_calcate_ratio

    Gₚ = Gᵖ(P, Z, sPON, gᶻ, p̃, kᶻ)
    
    return ((fᵈ * (1 - aᶻ) * (Gₚ + Gᵈ(P, Z, sPON, gᶻ, p̃, kᶻ))
            + fᵈ * mᵖ * P^2
            + (1 - fᶻ) * mᶻ * Z^2
            + (Gₚ * (1 - η) + mᵖ * P^2) * ρᶜᵃᶜᵒ³) * Rᵖ
            - μᵈᵈ * bPOC)
end

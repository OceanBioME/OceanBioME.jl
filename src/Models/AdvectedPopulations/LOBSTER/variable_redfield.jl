
@inline function (bgc::LOBSTER)(::Val{:sPOC}, x, y, z, t, NO₃, NH₄, P, Z, sPON, bPON, DOM, sPOC, dPOC, DOC, PAR)
    αᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    fᵈ = bgc.fast_sinking_mortality_fraction # really dumb definitions
    fᶻ = bgc.slow_sinking_mortality_fraction
    mᵖ = bgc.phytoplankton_mortality
    μᵈ = bgc.small_detritus_remineralisation_rate
    Rᵖ = bgc.phytoplankton_redfield

    return ((1 - fᵈ) * (1 - αᶻ) * (Gᵖ(P, Z, sPON, gᶻ, p̃, kᶻ) + Gᵈ(P, Z, sPON, gᶻ, p̃, kᶻ)) * Rᵖ
            + (1 - fᵈ) * mᵖ * P^2 * Rᵖ
            - Gᵈ(P, Z, sPON, gᶻ, p̃, kᶻ) * Rᵖ
            + fᶻ * mᶻ * Z^2 * Rᵖ
            - μᵈ * sPOC)
end

@inline function (bgc::LOBSTER)(::Val{:bPOC}, x, y, z, t, NO₃, NH₄, P, Z, sPON, bPON, DOM, sPOC, dPOC, DOC, PAR)
    αᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    fᵈ = bgc.fast_sinking_mortality_fraction # really dumb definitions
    fᶻ = bgc.slow_sinking_mortality_fraction
    mᵖ = bgc.phytoplankton_mortality
    μᵈᵈ = bgc.large_detritus_remineralisation_rate
    Rᵖ = bgc.phytoplankton_redfield

    return (fᵈ * (1 - αᶻ) * (Gᵖ(P, Z, sPON, gᶻ, p̃, kᶻ) + Gᵈ(P, Z, sPON, gᶻ, p̃, kᶻ)) * Rᵖ
            + fᵈ * mᵖ * P^2 * Rᵖ
            + (1 - fᶻ) * mᶻ * Z^2  * Rᵖ
            - μᵈᵈ * bPOC)
end
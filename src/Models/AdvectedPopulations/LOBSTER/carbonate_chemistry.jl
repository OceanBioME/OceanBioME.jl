@inline function (bgc::LOBSTER)(::Val{:DIC}, x, y, z, t, NO₃, NH₄, Pᶜ, Zᶜ, sPOC, bPOC, DOC, DIC, ALK, PAR)
    αᵖ = bgc.ammonia_fraction_of_exudate
    γ = bgc.phytoplankton_exudation_fraction
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    kₙₕ₄ = bgc.ammonia_half_saturation
    αᶻ = bgc.ammonia_fraction_of_excriment
    αᵈ = bgc.ammonia_fraction_of_detritus
    μᵈ = bgc.small_detritus_remineralisation_rate
    μᵈᵈ = bgc.large_detritus_remineralisation_rate
    μᵈᵒᵐ = bgc.disolved_organic_breakdown_rate
    μᶻ = bgc.zooplankton_excretion_rate
    ρᶜᵃᶜᵒ³ = bgc.organic_carbon_calcate_ratio

    return (μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * (αᵖ * γ - (1 + ρᶜᵃᶜᵒ³)) * Pᶜ
            + αᶻ * μᶻ * Zᶜ
            + αᵈ * μᵈ * sPOC
            + αᵈ * μᵈᵈ * bPOC
            + μᵈᵒᵐ * DOC)
end


@inline function (bgc::LOBSTER)(::Val{:ALK}, x, y, z, t, NO₃, NH₄, Pᶜ, Zᶜ, sPOC, bPOC, DOC, DIC, ALK, PAR)
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    kₙₕ₄ = bgc.ammonia_half_saturation
    ρᶜᵃᶜᵒ³ = bgc.organic_carbon_calcate_ratio

    return (μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) * Pᶜ
            - 2 * ρᶜᵃᶜᵒ³ * μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * Pᶜ)
end
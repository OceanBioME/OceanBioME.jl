
@inline function (bgc::LOBSTER)(::Val{:DIC}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR)
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
    Rᵈₚ = bgc.phytoplankton_redfield
    Rᵈₒ = bgc.disolved_organic_redfield
    ρᶜᵃᶜᵒ³ = bgc.organic_carbon_calcate_ratio

    return (μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * Rᵈₚ * (αᵖ * γ - (1 + ρᶜᵃᶜᵒ³)) * P
            + αᶻ * μᶻ * Rᵈₚ * Z
            + αᵈ * μᵈ * Dᶜ
            + αᵈ * μᵈᵈ * DDᶜ
            + μᵈᵒᵐ * DOM * Rᵈₒ)
end


@inline function (bgc::LOBSTER)(::Val{:ALK}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR)
    αᵖ = bgc.ammonia_fraction_of_exudate
    γ = bgc.phytoplankton_exudation_fraction
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    kₙₕ₄ = bgc.ammonia_half_saturation
    Rᵈₚ = bgc.phytoplankton_redfield
    ρᶜᵃᶜᵒ³ = bgc.organic_carbon_calcate_ratio

    return (μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) * P
            - 2 * ρᶜᵃᶜᵒ³ * μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * Rᵈₚ * P)
end
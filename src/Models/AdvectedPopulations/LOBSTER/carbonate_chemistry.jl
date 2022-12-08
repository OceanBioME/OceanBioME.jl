@inline function (bgc::LOBSTER)(::Val{:DIC}, x, y, z, t, NO₃, NH₄, P, Z, sPOM, bPOM, DOM, DIC, Alk, PAR)
    αᵖ = bgc.ammonia_fraction_of_exudate
    γ = bgc.phytoplankton_exudation_fraction
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    kₙₕ₄ = bgc.ammonia_half_saturation
    αᶻ = bgc.ammonia_fraction_of_excriment
    αᵈ = bgc.ammonia_fraction_of_detritus
    μˢᵖᵒᶜ = bgc.small_detritus_remineralisation_rate
    μᵇᵖᵒᶜ = bgc.large_detritus_remineralisation_rate
    μᵈᵒᵐ = bgc.disolved_organic_breakdown_rate
    μᶻ = bgc.zooplankton_excretion_rate
    ρᶜᵃᶜᵒ³ = bgc.organic_carbon_calcate_ratio
    Rdᵖ = bgc.phytoplankton_redfield
    Rdᵒ = bgc. organic_redfield
    η = bgc.zooplankton_calcite_dissolution
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation

    return (- μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * (1 + ρᶜᵃᶜᵒ³ * (1 - γ)) * P * Rdᵖ
            + αᵖ * γ * μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * P * Rdᵖ
            + αᶻ * μᶻ * Z * Rdᵖ
            + αᵈ * μˢᵖᵒᶜ * sPOM * Rdᵒ
            + αᵈ * μᵇᵖᵒᶜ * bPOM * Rdᵒ
            + μᵈᵒᵐ * DOM * Rdᵒ
            + Gᵖ(P, Z, sPOM, gᶻ, p̃, kᶻ) * (1 - η) * Rdᵖ * ρᶜᵃᶜᵒ³)
end


@inline function (bgc::LOBSTER)(::Val{:Alk}, x, y, z, t, NO₃, NH₄, P, Z, sPOM, bPOM, DOM, DIC, Alk, PAR)
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    kₙₕ₄ = bgc.ammonia_half_saturation
    ρᶜᵃᶜᵒ³ = bgc.organic_carbon_calcate_ratio
    Rd = bgc.phytoplankton_redfield

    return (μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) * P
            - 2 * ρᶜᵃᶜᵒ³ * μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * P * Rd)
end

# grazing
@inline p(P, D, p̃) = p̃ * P / (p̃ * P + (1 - p̃) * D + eps(0.0))
@inline Gᵈ(P, Z, D, gᶻ, p̃, kᶻ) = gᶻ * (1 - p(P, D, p̃)) * D * Z / (kᶻ + P * p(P, D, p̃) + (1 - p(P, D, p̃)) * D)
@inline Gᵖ(P, Z, D, gᶻ, p̃, kᶻ) = gᶻ * p(P, D, p̃) * P * Z / (kᶻ + P * p(P, D, p̃) + (1 - p(P, D, p̃)) * D)

# Limiting equations
@inline Lₚₐᵣ(PAR, kₚₐᵣ) = 1 - exp(-PAR/kₚₐᵣ)
@inline Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) = NO₃*exp(-ψ*NH₄)/(NO₃+kₙₒ₃)
@inline Lₙₕ₄(NH₄, kₙₕ₄) = max(0.0, NH₄/(NH₄+kₙₕ₄)) 

# Nutrients
@inline function (bgc::LOBSTER)(::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR)
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    μₙ = bgc.nitrifcaiton_rate

    return μₙ*NH₄ - μₚ*Lₚₐᵣ(PAR, kₚₐᵣ)*Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃)*P
end

@inline function (bgc::LOBSTER)(::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR)
    αᵖ = bgc.ammonia_fraction_of_exudate
    γ = bgc.phytoplankton_exudation_fraction
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    kₙₕ₄ = bgc.ammonia_half_saturation
    μₙ = bgc.nitrifcaiton_rate
    αᶻ = bgc.ammonia_fraction_of_excriment
    αᵈ = bgc.ammonia_fraction_of_detritus
    μᵈ = bgc.small_detritus_remineralisation_rate
    μᵈᵈ = bgc.large_detritus_remineralisation_rate
    μᵈᵒᵐ = bgc.disolved_organic_breakdown_rate
    μᶻ = bgc.zooplankton_excretion_rate

    return (αᵖ * γ * μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * P 
            - μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * Lₙₕ₄(NH₄, kₙₕ₄) * P
            - μₙ * NH₄
            + αᶻ * μᶻ * Z
            + αᵈ * μᵈ * D
            + αᵈ * μᵈᵈ * DD
            + μᵈᵒᵐ * DOM)
end

@inline function (bgc::LOBSTER)(::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR)
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

    return ((1 - αᵖ) * γ * μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * P 
            + (1 - αᶻ) * μᶻ * Z
            + (1 - αᵈ) * μᵈ * D
            + (1 - αᵈ) * μᵈᵈ * DD
            - μᵈᵒᵐ * DOM)
end

# Planktons
@inline function (bgc::LOBSTER)(::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR)
    γ = bgc.phytoplankton_exudation_fraction
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    kₙₕ₄ = bgc.ammonia_half_saturation
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᵖ = bgc.phytoplankon_mortality

    return ((1 - γ) * μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * (Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) + Lₙₕ₄(NH₄, kₙₕ₄)) * P 
            - Gᵖ(P, Z, D, gᶻ, p̃, kᶻ)
            - mᵖ * P^2)
end

@inline function (bgc::LOBSTER)(::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR)
    αᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    μᶻ = bgc.zooplankton_excretion_rate

    return (αᶻ * (Gᵖ(P, Z, D, gᶻ, p̃, kᶻ) + Gᵈ(P, Z, D, gᶻ, p̃, kᶻ))
            - mᶻ*Z^2
            - μᶻ*Z)
end

# Detritus

@inline function (bgc::LOBSTER)(::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR)
    αᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    fᵈ = bgc.fast_sinking_mortality_fraction # really dumb definitions
    fᶻ = bgc.slow_sinking_mortality_fraction
    mᵖ = bgc.phytoplankon_mortality
    μᵈ = bgc.small_detritus_remineralisation_rate

    return ((1 - fᵈ) * (1 - αᶻ) * (Gᵖ(P, Z, D, gᶻ, p̃, kᶻ) + Gᵈ(P, Z, D, gᶻ, p̃, kᶻ))
            + (1 - fᵈ) * mᵖ * P^2
            - Gᵈ(P, Z, D, gᶻ, p̃, kᶻ)
            + fᶻ * mᶻ * Z^2 
            - μᵈ * D)
end

@inline function (bgc::LOBSTER)(::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR)
    αᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    fᵈ = bgc.fast_sinking_mortality_fraction # really dumb definitions
    fᶻ = bgc.slow_sinking_mortality_fraction
    mᵖ = bgc.phytoplankon_mortality
    μᵈᵈ = bgc.large_detritus_remineralisation_rate

    return (fᵈ * (1 - αᶻ) * (Gᵖ(P, Z, D, gᶻ, p̃, kᶻ) + Gᵈ(P, Z, D, gᶻ, p̃, kᶻ))
            + fᵈ * mᵖ * P^2
            + (1 - fᶻ) * mᶻ * Z^2 
            - μᵈᵈ * DD)
end
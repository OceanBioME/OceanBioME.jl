using Oceananigans.Biogeochemistry: AbstractBiogeochemistry, maybe_velocity_fields
using Oceananigans.Units
using Oceananigans.Advection: CenteredSecondOrder

import Oceananigans.Biogeochemistry:
       required_biogeochemical_tracers,
       required_biogeochemical_auxiliary_fields,
       biogeochemical_drift_velocity,
       biogeochemical_advection_scheme


struct LOBSTER{FT, B, W, A} <: AbstractBiogeochemistry
    phytoplankton_preference :: FT
    maximum_grazing_rate :: FT
    grazing_half_saturation :: FT
    light_half_saturation :: FT
    nitrate_ammonia_inhibition :: FT
    nitrate_half_saturation :: FT
    ammonia_half_saturation :: FT
    maximum_phytoplankton_growthrate :: FT
    zooplankton_assimilation_fraction :: FT
    zooplankton_mortality :: FT
    zooplankton_excretion_rate :: FT
    phytoplankon_mortality :: FT
    small_detritus_remineralisation_rate :: FT
    large_detritus_remineralisation_rate :: FT
    phytoplankton_exudation_fraction :: FT
    nitrifcaiton_rate :: FT
    ammonia_fraction_of_exudate :: FT
    ammonia_fraction_of_excriment :: FT
    ammonia_fraction_of_detritus :: FT
    phytoplankton_redfield :: FT
    disolved_organic_redfield :: FT
    phytoplankton_chlorophyll_ratio :: FT
    organic_carbon_calcate_ratio :: FT
    respiraiton_oxygen_nitrogen_ratio :: FT
    nitrifcation_oxygen_nitrogen_ratio :: FT
    slow_sinking_mortality_fraction :: FT
    fast_sinking_mortality_fraction :: FT
    disolved_organic_breakdown_rate :: FT

    optionals :: B

    sinking_velocities :: W
    advection_schemes :: A

    function LOBSTER(;phytoplankton_preference::FT = 0.5,
                      maximum_grazing_rate::FT = 9.26e-6, # 1/s
                      grazing_half_saturation::FT = 1.0, # mmol N/m³
                      light_half_saturation::FT = 33.0, # W/m² (?)
                      nitrate_ammonia_inhibition::FT = 3.0,
                      nitrate_half_saturation::FT = 0.7, # mmol N/m³
                      ammonia_half_saturation::FT = 0.001, # mmol N/m³
                      maximum_phytoplankton_growthrate::FT = 1.21e-5, # 1/s
                      zooplankton_assimilation_fraction::FT = 0.7,
                      zooplankton_mortality::FT = 2.31e-6, # 1/s/mmol N/m³
                      zooplankton_excretion_rate::FT = 5.8e-7, # 1/s
                      phytoplankon_mortality::FT = 5.8e-7, # 1/s
                      small_detritus_remineralisation_rate::FT = 5.88e-7, # 1/s
                      large_detritus_remineralisation_rate::FT = 5.88e-7, # 1/s
                      phytoplankton_exudation_fraction::FT = 0.05,
                      nitrifcaiton_rate::FT = 5.8e-7, # 1/s
                      ammonia_fraction_of_exudate::FT = 0.75, 
                      ammonia_fraction_of_excriment::FT = 0.5,
                      ammonia_fraction_of_detritus::FT = 0.0,
                      phytoplankton_redfield::FT = 6.56, # mol C/mol N
                      disolved_organic_redfield::FT = 6.56, # mol C/mol N
                      phytoplankton_chlorophyll_ratio::FT = 1.31, # mol Chl/mol N
                      organic_carbon_calcate_ratio::FT = 0.1, # mol CaCO₃/mol N
                      respiraiton_oxygen_nitrogen_ratio::FT = 10.75, # mol O/molN
                      nitrifcation_oxygen_nitrogen_ratio::FT = 2.0, # mol O/molN
                      slow_sinking_mortality_fraction::FT = 0.5, 
                      fast_sinking_mortality_fraction::FT = 0.5,
                      disolved_organic_breakdown_rate::FT = 3.86e-7, # 1/s

                      carbonates = true,
                      oxygen = false,
                
                      sinking_velocities = (D = (0.0, 0.0, -3.47e-5), DD = (0.0, 0.0, -200/day)),
                      advection_schemes::A = NamedTuple{keys(sinking_velocities)}(repeat([CenteredSecondOrder()], 
                                                                                    length(sinking_velocities)))) where {FT, A}

        sinking_velocities = maybe_velocity_fields(sinking_velocities)
        W = typeof(sinking_velocities)
        optionals = Val((carbonates, oxygen))
        B = typeof(optionals)

        return new{FT, B, W, A}(phytoplankton_preference,
                                maximum_grazing_rate,
                                grazing_half_saturation,
                                light_half_saturation,
                                nitrate_ammonia_inhibition,
                                nitrate_half_saturation,
                                ammonia_half_saturation,
                                maximum_phytoplankton_growthrate,
                                zooplankton_assimilation_fraction,
                                zooplankton_mortality,
                                zooplankton_excretion_rate,
                                phytoplankon_mortality,
                                small_detritus_remineralisation_rate,
                                large_detritus_remineralisation_rate,
                                phytoplankton_exudation_fraction,
                                nitrifcaiton_rate,
                                ammonia_fraction_of_exudate,
                                ammonia_fraction_of_excriment,
                                ammonia_fraction_of_detritus,
                                phytoplankton_redfield,
                                disolved_organic_redfield,
                                phytoplankton_chlorophyll_ratio,
                                organic_carbon_calcate_ratio,
                                respiraiton_oxygen_nitrogen_ratio,
                                nitrifcation_oxygen_nitrogen_ratio,
                                slow_sinking_mortality_fraction,
                                fast_sinking_mortality_fraction,
                                disolved_organic_breakdown_rate,

                                optionals,
                            
                                sinking_velocities,
                                advection_schemes)
    end
end

@inline required_biogeochemical_tracers(::LOBSTER{<:Val{(false, false)}}) = (:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM)
@inline required_biogeochemical_tracers(::LOBSTER{<:Val{(true, false)}}) = (:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM, :DIC, :ALK)
@inline required_biogeochemical_tracers(::LOBSTER{<:Val{(false, true)}}) = (:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM, :OXY)
@inline required_biogeochemical_tracers(::LOBSTER{<:Val{(true, true)}}) = (:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM, :DIC, :ALK, :OXY)

@inline required_biogeochemical_auxiliary_fields(::LOBSTER) = (:PAR, )

const small_detritus = Union{Val{:D}, Val{:Dᶜ}}
const large_detritus = Union{Val{:DD}, Val{:DDᶜ}}

@inline biogeochemical_drift_velocity(bgc::LOBSTER, ::small_detritus) = bgc.sinking_velocity.D
@inline biogeochemical_drift_velocity(bgc::LOBSTER, ::large_detritus) = bgc.sinking_velocity.DD

@inline biogeochemical_advection_scheme(bgc::LOBSTER, ::small_detritus) = bgc.advection_scheme.D
@inline biogeochemical_advection_scheme(bgc::LOBSTER, ::large_detritus) = bgc.advection_scheme.DD

# grazing
@inline p(P, D, p̃) = p̃ * P / (p̃ * P + (1 - p̃) * D + eps(0.0))
@inline Gᵈ(P, Z, D, gᶻ, p̃, kᶻ) = gᶻ * (1 - p(P, D, p̃)) * D * Z / (kᶻ + P * p(P, D, p̃) + (1 - p(P, D, p̃)) * D)
@inline Gᵖ(P, Z, D, gᶻ, p̃, kᶻ) = gᶻ * p(P, D, p̃) * P * Z / (kᶻ + P * p(P, D, p̃) + (1 - p(P, D, p̃)) * D)

# Limiting equations
@inline Lₚₐᵣ(PAR, kₚₐᵣ) = 1 - exp(-PAR/kₚₐᵣ)
@inline Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) = NO₃*exp(-ψ*NH₄)/(NO₃+kₙₒ₃)
@inline Lₙₕ₄(NH₄, kₙₕ₄) = max(0.0, NH₄/(NH₄+kₙₕ₄)) 

# Nutrients
@inline function (bgc::LOBSTER)(::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    μₙ = bgc.nitrifcaiton_rate

    return μₙ*NH₄ - μₚ*Lₚₐᵣ(PAR, kₚₐᵣ)*Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃)*P
end

@inline function (bgc::LOBSTER)(::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
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


@inline function (bgc::LOBSTER)(::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
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
@inline function (bgc::LOBSTER)(::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
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

@inline function (bgc::LOBSTER)(::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
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

@inline function (bgc::LOBSTER)(::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
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

@inline function (bgc::LOBSTER)(::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
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

# Detritus carbon - these could be optional but could be complicated to add another option set with two optional variables
# as we won't be able to unequaly define forcing functions purely on the number of argumemnts
@inline function (bgc::LOBSTER)(::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
    αᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    fᵈ = bgc.fast_sinking_mortality_fraction # really dumb definitions
    fᶻ = bgc.slow_sinking_mortality_fraction
    mᵖ = bgc.phytoplankon_mortality
    μᵈ = bgc.small_detritus_remineralisation_rate
    Rᵈ = bgc.phytoplankton_redfield

    return (((1 - fᵈ) * (1 - αᶻ) * (Gᵖ(P, Z, D, gᶻ, p̃, kᶻ) + Gᵈ(P, Z, D, gᶻ, p̃, kᶻ))
            + (1 - fᵈ) * mᵖ * P^2
            - Gᵈ(P, Z, D, gᶻ, p̃, kᶻ)
            + fᶻ * mᶻ * Z^2) * Rᵈ
            - μᵈ * Dᶜ)
end

@inline function (bgc::LOBSTER)(::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
    αᶻ = bgc.zooplankton_assimilation_fraction
    gᶻ = bgc.maximum_grazing_rate
    p̃ = bgc.phytoplankton_preference
    kᶻ = bgc.grazing_half_saturation
    mᶻ = bgc.zooplankton_mortality
    fᵈ = bgc.fast_sinking_mortality_fraction # really dumb definitions
    fᶻ = bgc.slow_sinking_mortality_fraction
    mᵖ = bgc.phytoplankon_mortality
    μᵈᵈ = bgc.large_detritus_remineralisation_rate
    Rᵈ = bgc.phytoplankton_redfield

    return ((fᵈ * (1 - αᶻ) * (Gᵖ(P, Z, D, gᶻ, p̃, kᶻ) + Gᵈ(P, Z, D, gᶻ, p̃, kᶻ))
            + fᵈ * mᵖ * P^2
            + (1 - fᶻ) * mᶻ * Z^2) * Rᵈ
            - μᵈᵈ * DDᶜ)
end


# Carbonate chemistry
@inline (bgc::LOBSTER)(tracer::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)

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

# Oxygen chemistry
# there is a typo in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2010JC006446 so I am not sure the first term is correct, but this makes sense
@inline (bgc::LOBSTER)(tracer::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)

@inline function (bgc::LOBSTER)(::Val{:OXY}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR)
    γ = bgc.phytoplankton_exudation_fraction
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    ROᵖ = bgc.respiraiton_oxygen_nitrogen_ratio
    ROⁿ = bgc.nitrifcation_oxygen_nitrogen_ratio
    μₙ = bgc.nitrifcaiton_rate

    return (μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) * ROᵖ * P
            - (ROᴾ - ROⁿ)*bgc(Val(:NH₄), x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
            - ROᴾ * μₙ * NH₄)
end

# Oxygen and Carbonate
@inline (bgc::LOBSTER)(tracer::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DIC}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:ALK}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:OXY}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR)

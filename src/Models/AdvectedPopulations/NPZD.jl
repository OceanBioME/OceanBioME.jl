"""
Nutrient-Phytoplankton-Zooplankton-Detritus model of [Kuhn2015](@cite)
http://dx.doi.org/10.1016/j.pocean.2015.07.004
"""
module NPZD

export NutrientPhytoplanktonZooplanktonDetritus

using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Oceananigans.Advection: div_Uc
using Oceananigans.Forcings: maybe_constant_field
using Oceananigans.Units: day
using Oceananigans.Advection: WENO

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

struct NutrientPhytoplanktonZooplanktonDetritus <: AbstractBiogeochemistry
    # phytoplankton
    α # initial photosynthetic slope (1/(W/m²)/s)
    μ₀ # maximum growth at T=0°C (1/s)
    kₙ # Half-saturation coefficient of nutrient uptake (mmol N/m³)
    lᵖⁿ # respiration rate at T=0\degre°C (1/s)
    lᵖᵈ # mortality rarte at T=0°C (1/s)
    u⃗ᵖ # sinking rate (m/s)

    # zooplankton
    gₘₐₓ # maximum grazing rate (1/s)
    kₚ # half saturation coefficient of grazing (mmol N/m³)
    β # assimilation efficiency ()
    lᶻⁿ # excretion rate at T=0°C (1/s)
    lᶻᵈ # mortality rarte at T=0°C (1/s)

    # detritus
    rᵈⁿ # remineralization rate (1/s)
    u⃗ᵈ # sinking rate (m/s)

    # other things
    adv_scheme
end

# default parameters
function NutrientPhytoplanktonZooplanktonDetritus(; advection_scheme = WENO())
   #uᵖ, vᵖ, wᵖ = maybe_constant_field.((0.0, 0.0, -0.2551/day))
    u, v, w = maybe_constant_field.((0.0, 0.0, 0.0))
    u⃗ᵖ = (;u, v, w)

    #uᵈ, vᵈ, wᵈ = maybe_constant_field.((0.0, 0.0, -2.7489/day))
    u, v, w = maybe_constant_field.((0.0, 0.0, 0.0))
    u⃗ᵈ = (;u, v, w)

    return NutrientPhytoplanktonZooplanktonDetritus(0.1953/day, 0.6989/day, 2.3868, 0.066/day, 0.0101/day, u⃗ᵖ, 2.1522/day, sqrt(0.5573), 0.9116, 0.0102/day, 0.3395/day, 0.1213/day, u⃗ᵈ, advection_scheme)
end

function NutrientPhytoplanktonZooplanktonDetritus(α, μ₀, kₙ, lᵖⁿ, lᵖᵈ, wᴾ, gₘₐₓ, kₚ, β, lᶻⁿ, lᶻᵈ, rᵈⁿ, wᴰ; advection_scheme = WENO())
    u, v, w = maybe_constant_field.((0.0, 0.0, -wᴾ))
    u⃗ᵖ = (;u, v, w)

    u, v, w = maybe_constant_field.((0.0, 0.0, -wᴰ))
    u⃗ᵈ = (;u, v, w)

    return NutrientPhytoplanktonZooplanktonDetritus(α, μ₀, kₙ, lᵖⁿ, lᵖᵈ, u⃗ᵖ, gₘₐₓ, kₚ, β, lᶻⁿ, lᶻᵈ, rᵈⁿ, u⃗ᵈ, advection_scheme)
end

required_biogeochemical_tracers(::NutrientPhytoplanktonZooplanktonDetritus) = (:N, :P, :Z, :D, :T)  #:PAR) # we probably want to check for required auxiliary_fields
#=
required_biogeochemical_auxiliary_fields(::AbstractBiogeochemistry) = ()
required_biogeochemical_auxiliary_fields(::NutrientPhytoplanktonZooplanktonDetritus) = (:PAR, )
=#

@inline get_fields(i, j, k, names, fields) = @inbounds ntuple(n->fields[names][n][i, j, k], length(names))

@inline nutrient_limitation(N, kₙ) = N/(kₙ + N)

@inline Q₁₀(T) = 1.88^(T/10) # T in °C

@inline light_limitation(PAR, T, α, μ₀) = α*PAR/sqrt((μ₀*Q₁₀(T))^2 + α^2*PAR^2)

@inline function (bgc::NutrientPhytoplanktonZooplanktonDetritus)(i, j, k, grid, ::Val{:N}, clock, fields)
    N, P, Z, D, PAR, T = get_fields(i, j, k, (:N, :P, :Z, :D, :PAR, :T), fields)

    phytoplankton_consumption = bgc.μ₀*nutrient_limitation(N, bgc.kₙ)*light_limitation(PAR, T, bgc.α, bgc.μ₀)*P
    phytoplankton_metabolic_loss = bgc.lᵖⁿ*Q₁₀(T)*P
    zooplankton_metabolic_loss = bgc.lᶻⁿ*Q₁₀(T)*Z
    remineralization = bgc.rᵈⁿ*D

    return phytoplankton_metabolic_loss + zooplankton_metabolic_loss + remineralization - phytoplankton_consumption
end

@inline function (bgc::NutrientPhytoplanktonZooplanktonDetritus)(i, j, k, grid, ::Val{:P}, clock, fields)
    N, P, Z, PAR, T = get_fields(i, j, k, (:N, :P, :Z, :PAR, :T), fields)

    growth = bgc.μ₀*nutrient_limitation(N, bgc.kₙ)*light_limitation(PAR, T, bgc.α, bgc.μ₀)*P
    grazing = bgc.gₘₐₓ*nutrient_limitation(P^2, bgc.kₚ^2)*Z
    metabolic_loss = bgc.lᵖⁿ*Q₁₀(T)*P
    mortality_loss = bgc.lᵖᵈ*Q₁₀(T)*P

    sinking = div_Uc(i, j, k, grid, bgc.adv_scheme, bgc.u⃗ᵖ, fields.P)

    return growth - grazing - metabolic_loss - mortality_loss + sinking # hopefully sinking is signed right if we pass sinking speed as a positive number
end

@inline function (bgc::NutrientPhytoplanktonZooplanktonDetritus)(i, j, k, grid, ::Val{:P}, clock, fields)
    P, Z, T = get_fields(i, j, k, (:P, :Z, :T), fields)

    grazing = bgc.β*bgc.gₘₐₓ*nutrient_limitation(P^2, bgc.kₚ^2)*Z
    metabolic_loss = bgc.lᶻⁿ*Q₁₀(T)*Z
    mortality_loss = bgc.lᶻᵈ*Q₁₀(T)*Z^2

    return grazing - metabolic_loss - mortality_loss 
end

@inline function (bgc::NutrientPhytoplanktonZooplanktonDetritus)(i, j, k, grid, ::Val{:D}, clock, fields)
    P, Z, D, T = get_fields(i, j, k, (:P, :Z, :D, :T), fields)

    phytoplankton_mortality_loss = bgc.lᵖᵈ*Q₁₀(T)*P
    zooplankton_assimilation_loss = (1 - bgc.β)*bgc.gₘₐₓ*nutrient_limitation(P^2, bgc.kₚ^2)*Z
    zooplankton_mortality_loss = bgc.lᶻⁿ*Q₁₀(T)*Z^2

    remineralization = bgc.rᵈⁿ*D

    sinking = div_Uc(i, j, k, grid, bgc.adv_scheme, bgc.u⃗ᵈ, fields.D)

    return phytoplankton_mortality_loss + zooplankton_assimilation_loss + zooplankton_mortality_loss - remineralization + sinking
end

end
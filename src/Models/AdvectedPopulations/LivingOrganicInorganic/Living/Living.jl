module Living

struct Phytoplankton{FT, NL, LL, MO}
       base_growth_rate :: FT
    nutrient_limitation :: NL
       light_limitation :: LL
              mortality :: MO
end

@inline function (P::Phytoplankton)(i, j, k, grid, val_name::Val{:P}, clock, fields, auxiliary_fields)
    μ = P.base_growth_rate

    P = @inbounds fields.P[i, j, k]

    nutrient_limitaiton = P.nutrient_limitaiton(i, j, k, grid, clock, fields, auxiliary_fields)
    ligh_limitaiton     = P.ligh_limitaiton(i, j, k, grid, clock, fields, auxiliary_fields)
    mortality           = P.mortality(i, j, k, grid, clock, fields, auxiliary_fields)

    return μ * nutrient_limitaiton * light_limitation - mortality
end

struct Zooplankton{FT, TL, FC, MO, HG}
         base_grazing_rate :: FT
    temperature_limitation :: TL
        food_concentration :: FC
                 mortality :: MO
end

@inline function (Z::Zooplankton)(i, j, k, grid, val_name::Val{:Z}, clock, fields, auxiliary_fields)
    g = Z.base_grazing_rate

    Z = @inbounds fields.Z[i, j, k]

    temperature_limitation = Z.temperature_limitation(i, j, k, grid, clock, fields, auxiliary_fields)
    food_concentration     = P.food_concentration(i, j, k, grid, clock, fields, auxiliary_fields)
    mortality              = P.mortality(i, j, k, grid, clock, fields, auxiliary_fields)

    return g * food_concentration * Z - mortality
end

@inline function (Z::Zooplankton)(i, j, k, grid, val_name::Val{:P}, clock, fields, auxiliary_fields)
    

struct QuadraticHigherGrazing{FT}
    specific_rate :: FT
end

@inline (HG::QuadraticHigherGrazing)(i, j, k, grid, val_name::Val{:Z}, grid, clock, fields, auxiliary_fields) =
    @inbounds -HG.specific_rate *fields.Z[i, j, k]^2


@inline function (bgc::NPZD)(::Val{:P}, x, y, z, t, N, P, Z, D, T, PAR)
    μ₀ = bgc.base_maximum_growth
    kₙ = bgc.nutrient_half_saturation
    α = bgc.initial_photosynthetic_slope
    gₘₐₓ = bgc.maximum_grazing_rate
    kₚ = bgc.grazing_half_saturation
    lᵖⁿ = bgc.base_respiration_rate
    lᵖᵈ = bgc.phyto_base_mortality_rate

    growth = μ₀ * Q₁₀(T) * nutrient_limitation(N, kₙ) * light_limitation(PAR, α, μ₀ * Q₁₀(T)) * P
    grazing = gₘₐₓ * nutrient_limitation(P ^ 2, kₚ ^ 2) * Z
    metabolic_loss = lᵖⁿ * Q₁₀(T) * P
    mortality_loss = lᵖᵈ * Q₁₀(T) * P

    return growth - grazing - metabolic_loss - mortality_loss
end

@inline function (bgc::NPZD)(::Val{:Z}, x, y, z, t, N, P, Z, D, T, PAR)
    gₘₐₓ = bgc.maximum_grazing_rate
    kₚ = bgc.grazing_half_saturation
    lᶻⁿ = bgc.base_excretion_rate
    lᶻᵈ = bgc.zoo_base_mortality_rate
    β = bgc.assimulation_efficiency

    grazing = β * gₘₐₓ * nutrient_limitation(P ^ 2, kₚ ^ 2) * Z
    metabolic_loss = lᶻⁿ * Q₁₀(T) * Z
    mortality_loss = lᶻᵈ * Q₁₀(T) * Z ^ 2

    return grazing - metabolic_loss - mortality_loss
end


@inline function (lobster::PHYTO_ZOO_LOBSTER)(i, j, k, grid, val_name::Val{:P}, clock, fields, auxiliary_fields)
    γ = lobster.biology.phytoplankton_exudation_fraction
    m = lobster.biology.phytoplankton_mortality_rate

    P = @inbounds fields.P[i, j, k]

    μP = phytoplankton_growth(lobster, i, j, k, fields, auxiliary_fields)

    Gp = grazing(lobster, i, j, k, val_name, fields, auxiliary_fields)

    return (1 - γ) * μP - Gp - m * P^2
end # done!

@inline function (lobster::PHYTO_ZOO_LOBSTER)(i, j, k, grid, ::Val{:Z}, clock, fields, auxiliary_fields)
    α = lobster.biology.zooplankton_assimilation_fraction
    m = lobster.biology.zooplankton_mortality_rate
    μ = lobster.biology.zooplankton_excretion_rate

    Z = @inbounds fields.Z[i, j, k]
    G = total_grazing(lobster, i, j, k, fields, auxiliary_fields)

    return α * G - m * Z^2 - μ * Z
end # done!


end # module
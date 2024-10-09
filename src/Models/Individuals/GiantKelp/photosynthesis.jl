using Roots: find_zero, Bisection

@kwdef struct ArrheniusPhotosynthesis{FT}
       rate_at_ref_temp_1 :: FT = 1.22e-3 * 24
       rate_at_ref_temp_2 :: FT = 1.3e-3 * 24
               ref_temp_1 :: FT = 285.0
               ref_temp_2 :: FT = 288.0
           arrhenius_temp :: FT = (1 / ref_temp_1 - 1 / ref_temp_2) ^ -1 * log(rate_at_ref_temp_2 / rate_at_ref_temp_1)
                 low_temp :: FT = 271.0
                high_temp :: FT = 296.0
      high_arrhenius_temp :: FT = 1414.87
       low_arrhenius_temp :: FT = 4547.89
               efficiency :: FT = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60)
    saturation_irradiance :: FT = 90 * day/ (10 ^ 6)
      exudation_parameter :: FT = 0.5 # maybe this belongs somewhere else
  excess_exudate_redfield :: FT = Inf
end

@inline function (photo::ArrheniusPhotosynthesis)(T, PAR)
    PAR *= day / (3.99e-10 * 545e12) # W / m² / s to einstein / m² / day

    Tk = T + 273.15

    P₁ = photo.rate_at_ref_temp_1

    Tₐ  = photo.arrhenius_temp
    Tₐₗ = photo.low_arrhenius_temp
    Tₐₕ = photo.high_arrhenius_temp

    Tₚ  = photo.ref_temp_1
    Tₚₗ = photo.ref_temp_1
    Tₚₕ = photo.high_temp

    α = photo.efficiency
    Iₛ = photo.saturation_irradiance

    maximum_photosynthesis = P₁ * exp(Tₐ / Tₚ - Tₐ / Tk) / (1 + exp(Tₐₗ / Tk - Tₐₗ / Tₚₗ) + exp(Tₐₕ / Tₚₕ - Tₐₕ / Tk))
 
    β = solve_for_light_inhibition(photo, maximum_photosynthesis)

    pₛ = α * Iₛ / log(1 + α / β)

    return pₛ * (1 - exp(- α * PAR / pₛ)) * exp(-β * PAR / pₛ) 
end

@inline solve_for_light_inhibition(photo, maximum_photosynthesis) =
    find_zero(light_inhibition_residual, (0, 0.1), Bisection(); 
              p = (; maximum_photosynthesis, photo.efficiency, photo.saturation_irradiance))

@inline function light_inhibition_residual(β, p)
    pₘ = p.maximum_photosynthesis
    α = p.efficiency
    Iₛ = p.saturation_irradiance

    return pₘ - α * Iₛ / log(1 + α / β) * (α / (α + β)) * (β / (α + β)) ^ (β / α)
end

@inline function excudated_carbon_fraction(photo, kelp, C)
    γ = photo.exudation_parameter
    Cₘ = kelp.minimum_carbon_reserve

    return 1 - exp(γ * (Cₘ - C))
end

@inline function nitrogen_exudate(photo, kelp, C, T, PAR)
    CN = photo.excess_exudate_redfield

    P = photo(T, PAR)
    e = excudated_carbon_fraction(photo, kelp, C)

    return P * e * 14 / 12 / CN
end
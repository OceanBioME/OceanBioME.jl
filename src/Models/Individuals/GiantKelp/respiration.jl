@kwdef struct BasalAndUptakeRespiration{FT}
     base_basal_respiration_rate :: FT = 5.57e-5 * 24
  base_activity_respiration_rate :: FT = 1.11e-4 * 24

                      ref_temp_1 :: FT = 285.0
                      ref_temp_2 :: FT = 290.0

              rate_at_ref_temp_1 :: FT = 2.785e-4 * 24
              rate_at_ref_temp_2 :: FT = 5.429e-4 * 24

                  arrhenius_temp :: FT = (1 / ref_temp_1 - 1 / ref_temp_2) ^ -1 * log(rate_at_ref_temp_2 / rate_at_ref_temp_1)
end

@inline function total_respiration(resp, kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)
    Rᵇ = resp.base_basal_respiration_rate
    Rᵃ = resp.base_activity_respiration_rate

    Tₐ = resp.arrhenius_temp
    T₁ = resp.ref_temp_1

    Tk = T + 273.15

    f = exp(Tₐ / T₁ - Tₐ / Tk)

    μₘ = kelp.growth.maximum_specific_growth_rate
    Jₘ = kelp.nitrogen_uptake.maximum_nitrate_uptake + kelp.nitrogen_uptake.maximum_ammonia_uptake

    J = kelp.nitrogen_uptake(kelp, N, NO₃, NH₄, u, v, w) +
        ammonia_replacement(kelp.nitrogen_uptake, kelp, N, NH₄, u, v, w, μ)

    # (basal + activity associated) * temp factor
    return f * (Rᵇ + Rᵃ * (μ / μₘ + J / Jₘ))
end

@inline function (resp::BasalAndUptakeRespiration)(kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)
    Cₘ = kelp.minimum_carbon_reserve

    Rₜ = total_respiration(resp, kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)

    return Rₜ * (1 - exp(-(C - Cₘ) / Cₘ))
end

@inline function low_carbon_respiration(resp, kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)
    Cₘ = kelp.minimum_carbon_reserve

    Rₜ = total_respiration(resp, kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)

    return exp(-(C - Cₘ) / Cₘ) * Rₜ # this respiration product necessarily has redfield ratio Cₛ / Nₛ vs the normal which is unclear
end
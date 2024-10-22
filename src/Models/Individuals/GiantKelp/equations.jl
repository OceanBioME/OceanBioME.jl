@inline function (kelp::GiantKelp)(::Val{:A}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    Cₛ = kelp.structural_carbon

    μ = kelp.growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    ν = kelp.mortality_loss(kelp, t, A, N, C, u, v, w, T)

    R = low_carbon_respiration(kelp.respiration, kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)

    return A * (μ - ν - R / Cₛ) / day
end

@inline function (kelp::GiantKelp)(::Val{:N}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    kₐ = kelp.structural_dry_weight_per_area

    J = kelp.nitrogen_uptake(kelp, N, NO₃, NH₄, u, v, w)

    e = nitrogen_exudate(kelp.photosynthesis, kelp, C, T, PAR)

    consumption = nitrogen_consumption(kelp.growth, kelp, t, A, N, C, T, NH₄, u, v, w)

    return (J - e - consumption) / kₐ / day 
end

@inline function (kelp::GiantKelp)(::Val{:C}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    kₐ = kelp.structural_dry_weight_per_area

    P = kelp.photosynthesis(T, PAR)

    μ = kelp.growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    R = kelp.respiration(kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)

    e = excudated_carbon_fraction(kelp.photosynthesis, kelp, C)

    consumption = carbon_consumption(kelp.growth, kelp, t, A, N, C, T, NH₄, u, v, w)

    return ((1 - e) * P - R - consumption) / kₐ / day
end


#:NO₃, :NH₄, :DIC, :O₂, :DOC, :DON, :bPOC, :bPON

@inline function (kelp::SugarKelp)(::Val{:NO₃}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    J = nitrate_uptake(kelp, N, NO₃, u, v, w)

    return J * A / (day * 14 * 0.001) # gN/dm^2/hr to mmol N/s
end

@inline function (kelp::SugarKelp)(::Val{:NH₄}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    J = ammonia_uptake(kelp, t, A, N, C, T, NH₄, u, v, w)
    
    return J * A / (day * 14 * 0.001) # gN/dm^2/hr to mmol N/s
end

@inline function (kelp::SugarKelp)(::Val{:DIC}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    P = photosynthesis(kelp, T, PAR)

    μ = growth(kelp, t, A, N, C, T, NH₄, u, v, w)

    R = respiration(kelp, t, A, N, C, T, NO₃, NH₄, u, v, w, μ)
    
    return (P - R) * A / (day * 12 * 0.001)  # gC/dm^2/hr to mmol C/s
end 

# I now know that this may not be correct..., it probably should vary by where the nitrogen comes from
@inline (kelp::SugarKelp)(::Val{:O₂}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR) =
    - kelp(Val(:DIC), t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)

@inline function (kelp::SugarKelp)(::Val{:DOC}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    P = photosynthesis(kelp, T, PAR)

    e = specific_carbon_exudate(kelp, C)
    
    return e * P * A / (day * 12 * 0.001)  # gC/dm^2/hr to mmol C/s
end 

@inline (kelp::SugarKelp)(::Val{:DON}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR) =
    kelp(Val(:DOC), t, A, N, C, u, v, w, T, NO₃, NH₄, PAR) / kelp.exudation_redfield_ratio

@inline function (kelp::SugarKelp)(::Val{:bPON}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    kₐ = kelp.structural_dry_weight_per_area
    Nₛ = kelp.structural_nitrogen

    ν = erosion(kelp, t, A, N, C, T)

    return ν * kₐ * A * (N + Nₛ) / (day * 14 * 0.001)  # gN/dm^2/hr to mmol N/s
end 

# this is why the particles made by kelp can be much more carbon rich
@inline function (kelp::SugarKelp)(::Val{:bPOC}, t, A, N, C, u, v, w, T, NO₃, NH₄, PAR)
    kₐ = kelp.structural_dry_weight_per_area
    Cₛ = kelp.structural_carbon

    ν = erosion(kelp, t, A, N, C, T)

    return ν * kₐ * A * (C + Cₛ) / (day * 14 * 0.001)  # gC/dm^2/hr to mmol C/s
end 

struct LinearMortality{FT}
    mortality_rate :: FT
end

@inline (mort::LinearMortality)(kelp, t, A, N, C, u, v, w, T) = mort.mortality_rate

@kwdef struct SmallAreaLimitedMortality{FT}
     exponent :: FT = 0.22
    base_rate :: FT = 10^-6
end

@inline function (mort::SmallAreaLimitedMortality)(kelp, t, A, N, C, u, v, w, T)
    ε = mort.exponent
    ν₀ = mort.base_rate

    return ν₀ * exp(ε * A) / (1 + ν₀ * (exp(ε * A) - 1))
end
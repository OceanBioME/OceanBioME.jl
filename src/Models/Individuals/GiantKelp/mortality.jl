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

@kwdef struct SizeRegulatedMortality{FT}
    normalisation :: FT = 1.0
    base_rate :: FT = (0.0391+0.0199)/2
end

@inline function (mort::SizeRegulatedMortality)(kelp, t, A, N, C, u, v, w, T)
    A₀ = mort.normalisation
    ν₀ = mort.base_rate

    return ν₀ * tanh(A / A₀)
end

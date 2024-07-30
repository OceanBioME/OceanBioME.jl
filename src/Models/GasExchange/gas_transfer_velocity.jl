using Oceananigans.Units

"""
    ScaledTransferVelocity

k(660) transfer velocity parameterisation, defaults from Ho et al., 2006,
but `coeff` is wind product specific
"""
@kwdef struct ScaledTransferVelocity{FT, SC} 
             coeff :: FT = 0.266 / hour / 100 # cm/hour to m/s
    schmidt_number :: SC
end

(k::ScaledTransferVelocity)(u₁₀, T) = k.coeff * u₁₀ ^ 2 * (k.schmidt_number(T) / 660)^(-1/2)

Adapt.adapt_structure(to, k::ScaledTransferVelocity) = ScaledTransferVelocity(adapt(to, k.coeff),
                                                                              adapt(to, k.schmidt_number))
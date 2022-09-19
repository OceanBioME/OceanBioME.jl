@inline K_Fe(T) = 10^(16.27 - 1565.7/max(T, 5))
@inline Lₜ(DOC) = max(0.09*(DOC+40)-3.0, 0.6)

@inline function Feᶠ(Fe, DOC, T)
    K=K_Fe(T)
    Δ = 1 + K*Lₜ(DOC) - K*Fe
    return (-Δ + sqrt(Δ^2+4*K*Fe))/(2*K)
end

@inline function Cgfe1(DOC, POC, Fe, T, sh, a₁, a₂, a₄, a₅)
    FeL = Fe - Feᶠ(Fe, DOC, T)
    Fe_Coll = 0.5*FeL
    return ((a₁*DOC + a₂*POC)*sh + a₄*POC + a₅*DOC)*Fe_Coll
end

@inline function Cgfe2(GOC, sh, a₃)
    FeL = Fe - Feᶠ(Fe, DOC, T)
    Fe_Coll = 0.5*FeL
    return a₃*GOC*sh*Fe_Coll
end

@inline function BactfeR(μₚ, Lₗᵢₘᴮᵃᶜᵗ, bFe, K_Feᴮ¹)
    return μₚ*Lₗᵢₘᴮᵃᶜᵗ*θₘₐₓᶠᵉᴮᵃᶜᵗ*L_mondo(bFe, K_Feᴮ¹)
end
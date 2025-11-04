"""
    alkalinity_residual(H, p)

Returns the difference between total alkalinity computed from `H`` (hydrogen ion
concentration), `DIC`, `borate`, `sulfate`, `phosphate`, `silicate`, and `fluoride` 
concentration and chemical equilibrium constants specified in `p`, and the specified 
total `Alk`alinity.

    TAlk = [HCO₃⁻] + 2[CO₃²⁻] + [B(OH)₄⁻] + [OH⁻] + [HPO₄²⁻] + 2[PO₄³⁻] + [SiO(OH)₃⁻] 
           + [NH₃] + [HS⁻] - [H⁺] - [HSO₄⁻] - [HF] - [H₃PO₄] + minor acids and bases

Concentrations diagnosed as specified in Dickson et. al best practice descried in 
`CarbonChemistry` docstring.

Note ammonia (NH₃) is not currently included.
"""

@inline alkalinity_residual(H, p) =
     (bicarbonate(H, p)
     + carbonate(H, p)
     + borate(H, p)
     + hydroxide(H, p)
     + hydrogen_phosphate(H, p)
     + phosphate(H, p)
     + silicate(H, p)
     + free_hydrogen(H, p)
     + hydrogen_suplfate(H, p)
     + hydrogen_fluoride(H, p)
     + phosphoric_acid(H, p)
     - p.Alk)

@inline ∂ₕ_alkalinity_residual(H, p) =
    (∂ₕ_bicarbonate(H, p)
     + ∂ₕ_carbonate(H, p)
     + ∂ₕ_borate(H, p)
     + ∂ₕ_hydroxide(H, p)
     + ∂ₕ_hydrogen_phosphate(H, p)
     + ∂ₕ_phosphate(H, p)
     + ∂ₕ_silicate(H, p)
     + ∂ₕ_free_hydrogen(H, p)
     + ∂ₕ_hydrogen_suplfate(H, p)
     + ∂ₕ_hydrogen_fluoride(H, p)
     + ∂ₕ_phosphoric_acid(H, p))

carbonate_denom(H, p) = H^2 + p.K1 * H + p.K1 * p.K2
phosphorus_denom(H, p) = H^3 + p.KP1 * H^2 + p.KP1 * p.KP2 * H + p.KP1 * p.KP2 * p.KP3
sulfate_denom(H, p) = 1 + p.sulfate / p.KS

bicarbonate(H, p) =  p.K1 * H * p.DIC / carbonate_denom(H, p)
carbonate(H, p) = 2 * p.DIC * p.K1 * p.K2 / carbonate_denom(H, p)
borate(H, p) = p.boron / (1 + H / p.KB)
hydroxide(H, p) = p.KW / H
hydrogen_phosphate(H, p) = p.phosphate * p.KP1 * p.KP2 * H / phosphorus_denom(H, p)
phosphate(H, p) = 2 * p.phosphate * p.KP1 * p.KP2 * p.KP3 / phosphorus_denom(H, p)
silicate(H, p) = p.silicate / (1 + H / p.KSi)
free_hydrogen(H, p) = - H / sulfate_denom(H, p)
hydrogen_suplfate(H, p) = - p.sulfate / (1 + p.KS / H * sulfate_denom(H, p))
hydrogen_fluoride(H, p) = -p.fluoride / (1 + p.KF / H)
phosphoric_acid(H, p) = -p.phosphate * H^3 / phosphorus_denom(H, p)

∂ₕ_carbonate_denom(H, p) = 2 * H + p.K1 
∂ₕ_phosphorus_denom(H, p) = 3 * H^2 + 2 * p.KP1 * H + p.KP1 * p.KP2
∂ₕ_sulfate_denom(H, p) = 0

∂ₕ_bicarbonate(H, p) = p.K1 * p.DIC * (p.K1*p.K2 - H^2) / carbonate_denom(H, p)^2
∂ₕ_carbonate(H, p) = -∂ₕ_carbonate_denom(H, p) * 2 * p.DIC * p.K1 * p.K2 / (carbonate_denom(H, p)^2)
∂ₕ_borate(H, p) = -p.boron / (1 + H / p.KB)^2 / p.KB
∂ₕ_hydroxide(H, p) = -p.KW / H^2
∂ₕ_hydrogen_phosphate(H, p) = p.phosphate * p.KP1 * p.KP2 / phosphorus_denom(H, p) - ∂ₕ_phosphorus_denom(H, p) * p.phosphate * p.KP1 * p.KP2 * H / phosphorus_denom(H, p)^2
∂ₕ_phosphate(H, p) = - ∂ₕ_phosphorus_denom(H, p) * 2 * p.phosphate * p.KP1 * p.KP2 * p.KP3 / phosphorus_denom(H, p)^2
∂ₕ_silicate(H, p) = -p.silicate / (1 + H / p.KSi)^2 / p.KSi
∂ₕ_free_hydrogen(H, p) = - 1 / sulfate_denom(H, p) 
∂ₕ_hydrogen_suplfate(H, p) = p.sulfate / (1 + p.KS / H * sulfate_denom(H, p))^2 * (p.KS / H * ∂ₕ_sulfate_denom(H, p) - p.KS / H^2 * sulfate_denom(H, p))
∂ₕ_hydrogen_fluoride(H, p) = -p.fluoride / (1 + p.KF / H) ^ 2 * p.KF / H^2
∂ₕ_phosphoric_acid(H, p) = -3 * p.phosphate * H^2 / phosphorus_denom(H, p) - ∂ₕ_phosphorus_denom(H, p) * p.phosphate * H^3 / phosphorus_denom(H, p)^2 

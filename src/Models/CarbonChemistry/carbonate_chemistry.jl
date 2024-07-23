using Roots
using OceanBioME.Models: teos10_polynomial_approximation

"""
    CarbonChemistry(; ionic_strength = IonicStrength(),
                      solubility = K0(),
                      carbonic_acid = (K1 = K1(), K2 = K2()),
                      boric_acid = KB(),
                      water = KW(),
                      sulfate = KS(; ionic_strength),
                      fluoride = KF(; ionic_strength),
                      phosphoric_acid = (KP1 = KP1(), KP2 = KP2(), KP3 = KP3()),
                      silicic_acid = KSi(; ionic_strength),
                      calcite_solubility = KSP_calcite(),
                      density_function = teos10_polynomial_approximation)

Carbon chemistry model capable of solving for sea water pCO₂ from DIC and 
total alkalinity or DIC and pH. 

Default form from Dickson, A.G., Sabine, C.L. and  Christian, J.R. (2007), 
Guide to Best Practices for Ocean CO 2 Measurements. PICES Special Publication 3, 191 pp.

See each parameters documentation for origional sources.

Example
=======

```jldoctest
julia> using OceanBioME

julia> carbon_chemistry = CarbonChemistry()
Carbon chemistry model which solves for pCO₂ and pH

julia> pCO₂ = carbon_chemistry(2000, 2000, 10, 35)
1325.07184935965

julia> pH = carbon_chemistry(2000, 2000, 10, 35; return_pH = true)
7.494548148653669

julia> pCO₂_higher_pH = carbon_chemistry(2000, NaN, 10, 35, 7.5)
1308.7134398483074

```
"""
@kwdef struct CarbonChemistry{P0, PC, PB, PS, PF, PP, PSi, PW, IS, PKS, PRho}
          ionic_strength :: IS   = IonicStrength()
              solubility :: P0   = K0()
           carbonic_acid :: PC   = (K1 = K1(), K2 = K2())
              boric_acid :: PB   = KB()
                   water :: PW   = KW()
                 sulfate :: PS   = KS(; ionic_strength)
                fluoride :: PF   = KF(; ionic_strength)
         phosphoric_acid :: PP   = (KP1 = KP1(), KP2 = KP2(), KP3 = KP3())
            silicic_acid :: PSi  = KSi(; ionic_strength)
      calcite_solubility :: PKS  = KSP_calcite()
        density_function :: PRho = teos10_polynomial_approximation
end

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
@inline function alkalinity_residual(H, p)
    carbonate_denom = H^2 + p.K1 * H + p.K1 * p.K2
    phosphorus_denom = H^3 + p.KP1 * H^2 + p.KP1 * p.KP2 * H + p.KP1 * p.KP2 * p.KP3
    sulfate_denom = 1 + p.sulfate / p.KS

    bicarbonate = p.K1 * H * p.DIC / carbonate_denom
    carbonate = 2 * p.DIC * p.K1 * p.K2 / carbonate_denom
    borate = p.boron / (1 + H / p.KB)
    hydroxide = p.KW / H
    hydrogen_phosphate = p.phosphate * p.KP1 * p.KP2 * H / phosphorus_denom
    phosphate = 2 * p.phosphate * p.KP1 * p.KP2 * p.KP3 / phosphorus_denom
    silicate = p.silicate / (1 + H / p.KSi)
    free_hydrogen = - H / sulfate_denom
    hydrogen_suplfate = - p.sulfate / (1 + p.KS / H / sulfate_denom)
    hydrogen_fluoride = -p.fluoride / (1 + p.KF / H)
    phosphoric_acid = -p.phosphate * H^3 / phosphorus_denom

    return (bicarbonate 
            + carbonate
            + borate
            + hydroxide
            + hydrogen_phosphate
            + phosphate
            + silicate
            + free_hydrogen
            + hydrogen_suplfate
            + hydrogen_fluoride
            + phosphoric_acid 
            - p.Alk)
end

"""
    (p::CarbonChemistry)(; DIC, T, S, Alk = 0, pH = nothing,
                           return_pH = false,
                           boron = 0.000232 / 10.811 * S / 1.80655,
                           sulfate = 0.14 / 96.06 * S / 1.80655,
                           fluoride = 0.000067 / 18.9984 * S / 1.80655,
                           silicate = 0,
                           phosphate = 0,
                           upper_pH_bound = 14,
                           lower_pH_bound = 0)

Calculates `pCO₂` in sea water with `DIC`, `Alk`alinity, `T`emperature, and `S`alinity
unless `pH` is specified, in which case intermediate computation of `pH` is skipped and
`pCO₂` is calculated from the `DIC`, `T`, `S` and `pH`.

Alternativly, `pH` is returned if `return_pH` is `true`.
"""
@inline function (p::CarbonChemistry)(; DIC, T, S, Alk = 0, pH = nothing,
                                        P = nothing,
                                        lon = 0,
                                        lat = 0,
                                        return_pH = false,
                                        boron = 0.000232 / 10.811 * S / 1.80655,
                                        sulfate = 0.14 / 96.06 * S / 1.80655,
                                        fluoride = 0.000067 / 18.9984 * S / 1.80655,
                                        silicate = 0,
                                        phosphate = 0,
                                        upper_pH_bound = 14,
                                        lower_pH_bound = 0)

    ρₒ = p.density_function(T, S, ifelse(isnothing(P), 0, P), lon, lat)

    # Centigrade to kelvin
    T += 273.15

    # mili-equivalents / m³ to equivalents / kg
    Alk *= 1e-3 / ρₒ

    # mmol / m³ to mol / kg
    DIC       *= 1e-3 / ρₒ
    phosphate *= 1e-3 / ρₒ
    silicate  *= 1e-3 / ρₒ
    
    # ionic strength
    Is = p.ionic_strength(S)

    # compute equilibrium constants
    K1  = p.carbonic_acid.K1(T, S; P)
    K2  = p.carbonic_acid.K2(T, S; P)
    KB  = p.boric_acid(T, S; P)
    KW  = p.water(T, S; P)
    KS  = p.sulfate(T, S, Is; P)
    KF  = p.fluoride(T, S, Is, KS; P)
    KP1 = p.phosphoric_acid.KP1(T, S; P)
    KP2 = p.phosphoric_acid.KP2(T, S; P)
    KP3 = p.phosphoric_acid.KP3(T, S; P)
    KSi = p.silicic_acid(T, S, Is)

    params = (; DIC, Alk, boron, sulfate, fluoride, silicate, phosphate,
                K1, K2, KB, KW, KS, KF, KP1, KP2, KP3, KSi)

    # solve equilibrium for hydrogen ion concentration
    H = solve_for_H(pH, params, upper_pH_bound, lower_pH_bound)

    # compute solubility equilibrium constant
    K0 = p.solubility(T, S)

    # compute pCO₂
    CO₂  = DIC * H ^ 2 / (H ^ 2 + K1 * H + K1 * K2) 
    pCO₂ = (CO₂ / K0) * 10 ^ 6 # μatm

    # compute pH
    pH = -log10(H)

    return ifelse(return_pH, pH, pCO₂) 
end

# solves `alkalinity_residual` for pH
solve_for_H(pH, args...) = 10.0 ^ - pH

solve_for_H(::Nothing, params, upper_pH_bound, lower_pH_bound) =
    find_zero(alkalinity_residual, (10.0 ^ - upper_pH_bound, 10.0 ^ - lower_pH_bound), Bisection(); atol = 1e-10, p = params)

# display
summary(::IO, ::CarbonChemistry) = string("Carbon chemistry model")
show(io::IO, ::CarbonChemistry) = print(io, "Carbon chemistry model which solves for pCO₂ and pH")
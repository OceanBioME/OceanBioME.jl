"""
    carbonate_concentration(cc::CarbonChemistry; DIC, T, S, Alk, ...)

Compute the carbonate ion concentration, `[CO₃²⁻]`, by solving the carbon chemistry
equilibrium for the hydrogen ion concentration and speciating `DIC`.

`DIC` is expected in mmol C/m³, `Alk` in meq/m³, `silicate` and `phosphate` in mmol/m³,
`T` in °C, `S` in PSU, and `P` in bar; `boron`, `sulfate`, and `fluoride` are in mol/kg.
The returned `[CO₃²⁻]` is in mmol/m³, matching `DIC`.

When `pH` is specified the intermediate solve is skipped and the free pH (i.e. -log[H⁺])
is expected.
"""
function carbonate_concentration(cc::CarbonChemistry;
                                 DIC::FT, T, S, Alk = zero(DIC), pH = nothing,
                                 P = nothing,
                                 lon = zero(DIC),
                                 lat = zero(DIC),
                                 boron = convert(typeof(DIC), 0.000232 / 10.811 * S / 1.80655),
                                 sulfate = convert(typeof(DIC), 0.14 / 96.06 * S / 1.80655),
                                 fluoride = convert(typeof(DIC), 0.000067 / 18.9984 * S / 1.80655),
                                 silicate = zero(DIC),
                                 phosphate = zero(DIC),
                                 initial_pH_guess = convert(typeof(DIC), 8)) where FT

    ρₒ = cc.density_function(T, S, ifelse(isnothing(P), zero(DIC), P), lon, lat)

    # Centigrade to kelvin
    T += convert(FT, 273.15)

    # mili-equivalents / m³ to equivalents / kg
    Alk *= convert(FT, 1e-3) / ρₒ

    # mmol / m³ to mol / kg
    DIC       *= convert(FT, 1e-3) / ρₒ
    phosphate *= convert(FT, 1e-3) / ρₒ
    silicate  *= convert(FT, 1e-3) / ρₒ
    
    # ionic strength
    Is = cc.ionic_strength(S)

    # compute equilibrium constants
    K1 = cc.carbonic_acid.K1(T, S; P)
    K2 = cc.carbonic_acid.K2(T, S; P)
    KB = cc.boric_acid(T, S; P)
    KW = cc.water(T, S; P)
    KS = cc.sulfate(T, S, Is; P)
    KF = cc.fluoride(T, S, Is, KS; P)
    KP1 = cc.phosphoric_acid.KP1(T, S; P)
    KP2 = cc.phosphoric_acid.KP2(T, S; P)
    KP3 = cc.phosphoric_acid.KP3(T, S; P)
    KSi = cc.silicic_acid(T, S, Is; P)

    params = (; DIC, Alk, boron, sulfate, fluoride, silicate, phosphate,
                K1, K2, KB, KW, KS, KF, KP1, KP2, KP3, KSi)

    # solve equilibrium for hydrogen ion concentration

    H = solve_for_H(pH, params, initial_pH_guess, cc.solver)

    # compute the carbonate ion concentration (mol / kg)
    denom1 = (H * (H + K1))
    denom2 = (one(DIC) + K1 * K2 / denom1)

    CO₃²⁻ = DIC * K1 * K2 / denom1 / denom2

    # mol / kg to mmol / m³
    return CO₃²⁻ * ρₒ * convert(FT, 1e3)
end

"""
    calcium_carbonate_saturation(cc::CarbonChemistry; DIC, T, S, Alk, ...)

Compute the calcium carbonate saturation state, `Ω = [Ca²⁺][CO₃²⁻] / KSP`, where `KSP` is
the solubility of the mineral phase set by `cc.calcium_carbonate_solubility` (calcite by
default).

Units follow [`carbonate_concentration`](@ref): `DIC` in mmol C/m³, `Alk` in meq/m³,
`silicate` and `phosphate` in mmol/m³, `T` in °C, `S` in PSU, and `P` in bar, while
`calcium_ion_concentration` (like `boron`, `sulfate`, and `fluoride`) is in mol/kg since
that is the basis `KSP` is defined on. `Ω` is dimensionless.
"""
function calcium_carbonate_saturation(cc::CarbonChemistry;
                            DIC::FT, T, S, Alk = zero(DIC), pH = nothing,
                            P = nothing,
                            lon = zero(DIC),
                            lat = zero(DIC),
                            boron = convert(typeof(DIC), 0.000232 / 10.811 * S / 1.80655),
                            sulfate = convert(typeof(DIC), 0.14 / 96.06 * S / 1.80655),
                            fluoride = convert(typeof(DIC), 0.000067 / 18.9984 * S / 1.80655),
                            calcium_ion_concentration = convert(typeof(DIC), 0.0103 * S / 35),
                            silicate = zero(DIC),
                            phosphate = zero(DIC),
                            initial_pH_guess = convert(typeof(DIC), 8)) where FT

    CO₃²⁻ = carbonate_concentration(cc;
                                    DIC, Alk, T, S, pH,
                                    P,
                                    lon,
                                    lat,
                                    boron,
                                    sulfate,
                                    fluoride,
                                    silicate,
                                    phosphate,
                                    initial_pH_guess)

    ρₒ = cc.density_function(T, S, ifelse(isnothing(P), zero(DIC), P), lon, lat)

    # KSP is defined on a mol / kg basis, so put [CO₃²⁻] back on one too
    CO₃²⁻ *= convert(FT, 1e-3) / ρₒ

    KSP = cc.calcium_carbonate_solubility(T + convert(FT, 273.15), S; P)

    return calcium_ion_concentration * CO₃²⁻ / KSP
end

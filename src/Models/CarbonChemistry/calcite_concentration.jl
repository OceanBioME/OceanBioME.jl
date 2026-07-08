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

    # compute the calcite concentration
    denom1 = (H * (H + K1))
    denom2 = (one(DIC) + K1 * K2 / denom1)

    return DIC * K1 * K2 / denom1 / denom2
end

function calcite_saturation(cc::CarbonChemistry;
                            DIC::FT, T, S, Alk = zero(DIC), pH = nothing,
                            P = nothing,
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
                                    boron,
                                    sulfate,
                                    fluoride,
                                    silicate,
                                    phosphate,
                                    initial_pH_guess)

    KSP = cc.calcite_solubility(T+convert(FT, 273.15), S; P)

    # not confident these all have the right units
    return calcium_ion_concentration * CO₃²⁻ / KSP
end

"""
    carbon_dioxide_and_carbonate(cc::CarbonChemistry; DIC, T, S, Alk, P, …)

Return the aqueous carbon-dioxide (`[CO₂(aq)] = [H₂CO₃]`) and carbonate-ion (`[CO₃²⁻]`)
concentrations, both in **mmol m⁻³** (i.e. the same units as `DIC`/`Alk`, and as MARBL's
`carbonate%H2CO3`/`CO3`), from a single pH solve. Mirrors [`carbonate_concentration`] but returns
CO₂ alongside CO₃ and converts back from the internal mol/kg to mmol m⁻³ (`× ρ × 10³`). Used by
`ExplicitCalcite` to precompute the `:CO₂`/`:CO₃` auxiliary fields the biology reads.

TODO: remove this redundancy and make units optional
"""
function carbon_dioxide_and_carbonate(cc::CarbonChemistry;
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

    # mili-equivalents / m³ to equivalents / kg, mmol / m³ to mol / kg
    Alk       *= convert(FT, 1e-3) / ρₒ
    DIC       *= convert(FT, 1e-3) / ρₒ
    phosphate *= convert(FT, 1e-3) / ρₒ
    silicate  *= convert(FT, 1e-3) / ρₒ

    Is = cc.ionic_strength(S)

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

    H = solve_for_H(pH, params, initial_pH_guess, cc.solver)

    # speciate DIC into CO₂(aq) and CO₃²⁻ (mol/kg), then back to mmol/m³
    denom = H ^ 2 + K1 * H + K1 * K2
    CO₂   = DIC * H ^ 2 / denom      # [CO₂(aq)] = [H₂CO₃]
    CO₃   = DIC * K1 * K2 / denom    # [CO₃²⁻]

    to_volume = ρₒ * convert(FT, 1e3)   # mol/kg → mmol/m³

    return CO₂ * to_volume, CO₃ * to_volume
end
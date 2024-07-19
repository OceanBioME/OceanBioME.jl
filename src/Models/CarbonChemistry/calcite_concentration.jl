function calcite_concentration(cc::CarbonChemistry, 
                               DIC, Alk, T, S, pH = nothing;
                               P = nothing,
                               boron = 0.000232 / 10.811 * S / 1.80655,
                               sulfate = 0.14 / 96.06 * S / 1.80655,
                               fluoride = 0.000067 / 18.9984 * S / 1.80655,
                               silicate = 0,
                               phosphate = 0,
                               upper_pH_bound = 14,
                               lower_pH_bound = 0)

    
    ρₒ = seawater_density(T, S)

    # Centigrade to kelvin
    T += 273.15

    # mili-equivalents / m³ to equivalents / kg
    Alk *= 1e-3 / ρₒ

    # mmol / m³ to mol / kg
    DIC       *= 1e-3 / ρₒ
    phosphate *= 1e-3 / ρₒ
    silicate  *= 1e-3 / ρₒ
    
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
    KSi = cc.silicic_acid(T, S, Is)

    params = (; DIC, Alk, boron, sulfate, fluoride, silicate, phosphate,
                K1, K2, KB, KW, KS, KF, KP1, KP2, KP3, KSi)

    # solve equilibrium for hydrogen ion concentration
    H = solve_for_H(pH, params, upper_pH_bound, lower_pH_bound)

    # compute the calcite concentration
    denom1 = (H * (H + K1))
    denom2 = (1 + K1 * K2 / denom1)

    return DIC * K1 * K2 / denom1 / denom2
end

function calcite_saturation(cc::CarbonChemistry, 
                            DIC, Alk, T, S, pH = nothing;
                            P = nothing,
                            boron = 0.000232 / 10.811 * S / 1.80655,
                            sulfate = 0.14 / 96.06 * S / 1.80655,
                            fluoride = 0.000067 / 18.9984 * S / 1.80655,
                            calcium_ion_concentration = 0.0103 * S / 35,
                            silicate = 0,
                            phosphate = 0,
                            upper_pH_bound = 14,
                            lower_pH_bound = 0)

    CO₃²⁻ = calcite_concentration(cc, 
                                  DIC, Alk, T, S, pH;
                                  P,
                                  boron,
                                  sulfate,
                                  fluoride,
                                  silicate,
                                  phosphate,
                                  upper_pH_bound,
                                  lower_pH_bound)

    KSP = cc.calcite_solubility(T, S; P)

    return calcium_ion_concentration * CO₃²⁻ / KSP
end
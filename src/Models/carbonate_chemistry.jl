##### 
##### Carbonate chemistry to determine pCO₂
##### As per OCMIP Phase 2 http://ocmip5.ipsl.jussieu.fr/OCMIP/phase2/simulations/Abiotic/Cchem/co2calc.f simplified as in PISCESv2
#####

using Roots

"""
    seawater_density(T, S)

Returns the density of seawater at `T` °C and `S` psu in kg/m³ as
per Chapter 5, Section 4.2 of Dickson, Sabine and Christian (2007, 
http://cdiac.ornl.gov/oceans/Handbook_2007.html).

With many thanks to Oscar Branson for his implementation in cbsyst
(https://github.com/oscarbranson/cbsyst).
"""
@inline function seawater_density(T, S)
    # convert temperature to IPTS-68
    T = (T + 0.0002) / 0.99975

    pSMOW = (999.842594
             + 6.793952e-2 * T
             - 9.095290e-3 * T^2
             + 1.001685e-4 * T^3
             - 1.120083e-6 * T^4
             + 6.536332e-9 * T^5)

    A = (8.24493e-1
         - 4.0899e-3 * T
         + 7.6438e-5 * T^2
         - 8.2467e-7 * T^3
         + 5.3875e-9 * T^4)

    B = -5.72466e-3 + 1.0227e-4 * T - 1.6546e-6 * T^2
    C = 4.8314e-4

    return pSMOW + A * S + B * S^1.5 + C * S^2
end

include("equilibrium_constants.jl")

@kwdef struct CarbonChemistry{P0, PC, PB, PS, PF, PP, PSi, PW, IS}
          ionic_strength :: IS  = IonicStrength()
              solubility :: P0  = K0()
           carbonic_acid :: PC  = (K1 = K1(), K2 = K2())
              boric_acid :: PB  = KB()
                   water :: PW  = KW()
                 sulfate :: PS  = KS(; ionic_strength)
                fluoride :: PF  = KF(; ionic_strength)
         phosphoric_acid :: PP  = (KP1 = KP1(), KP2 = KP2(), KP3 = KP3())
            silicic_acid :: PSi = KSi(; ionic_strength)
end

@inline function solve_for_H(H, p)
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

@inline function (p::CarbonChemistry)(DIC, Alk, T, S;
                                      pH = nothing,
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

    # mili-equivilants / m³ to equivilants / kg
    Alk *= 1e-3 / ρₒ

    # mmol / m³ to mol / kg
    DIC *= 1e-3 / ρₒ
    phosphate *= 1e-3 / ρₒ
    silicate *= 1e-3 / ρₒ
    
    # ionic strength
    Is = p.ionic_strength(S)

    K1 = p.carbonic_acid.K1(T, S)
    K2 = p.carbonic_acid.K2(T, S)
    KB = p.boric_acid(T, S)
    KW = p.water(T, S)
    KS = p.sulfate(T, S, Is)
    KF = p.fluoride(T, S, Is, KS)
    KP1 = p.phosphoric_acid.KP1(T, S)
    KP2 = p.phosphoric_acid.KP2(T, S)
    KP3 = p.phosphoric_acid.KP3(T, S)
    KSi = p.silicic_acid(T, S, Is)

    params = (; DIC, Alk, boron, sulfate, fluoride, silicate, phosphate,
                K1, K2, KB, KW, KS, KF, KP1, KP2, KP3, KSi)

    H = compute_H(pH, params, upper_pH_bound, lower_pH_bound)

    FF = p.solubility(T, S)

    CO₂ = DIC * H ^ 2 / (H ^ 2 + K1 * H + K1 * K2) 
    pCO₂ = (CO₂ / FF) * 10 ^ 6

    return -log10(H)#pCO₂ # μatm
end

compute_H(pH, args...) = 10.0 ^ - pH

compute_H(pH::Nothing, params, upper_pH_bound, lower_pH_bound) =
    find_zero(solve_for_H, (10.0 ^ - upper_pH_bound, 10.0 ^ - lower_pH_bound), Bisection(); atol = 1e-10, p = params)
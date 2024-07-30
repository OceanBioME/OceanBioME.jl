"""
    CarbonDioxideConcentration(; carbon_chemistry, 
                                 first_virial_coefficient = PolynomialVirialCoefficientForCarbonDioxide(),
                                 cross_viral_coefficient = CrossVirialCoefficientForCarbonDioxide(),
                                 air_pressue = 1 # atm)

Converts fCO₂ to partial pressure as per Dickson, A.G., Sabine, C.L. and  Christian, J.R. (2007), 
Guide to Best Practices for Ocean CO 2 Measurements. PICES Special Publication 3, 191 pp.
"""
@kwdef struct CarbonDioxideConcentration{CC<:CarbonChemistry, FV, CV, AP} <: Function
            carbon_chemistry :: CC 
    first_virial_coefficient :: FV = PolynomialVirialCoefficientForCarbonDioxide()
    cross_virial_coefficient :: CV = CrossVirialCoefficientForCarbonDioxide()
                air_pressure :: AP = 1 # atm
end

field_dependencies(::CarbonDioxideConcentration) = (:T, :S, :DIC, :Alk)
optional_fields(::CarbonDioxideConcentration) = (:silicate, :phosphate)

@inline function (cc::CarbonDioxideConcentration)(x, y, t, T, S, DIC, Alk, args...)
    fCO₂ = call_carbon_chemistry(cc.carbon_chemistry, DIC, Alk, T, S, args...)
    
    P = surface_value(cc.air_pressure, x, y, t) * ATM
    Tk = T + 273.15

    B = cc.first_virial_coefficient(Tk)
    δ = cc.cross_virial_coefficient(Tk)

    xCO₂ = fCO₂ / P
    
    # Experimentally this converged xCO₂ to machine precision
    for n = 1:3
        φ = P * exp((B + 2 * (1 - xCO₂)^2 * δ) * P / (GAS_CONSTANT * Tk))
        xCO₂ = fCO₂ / φ
    end

    return xCO₂ * P # mol/mol to ppm
end

@inline call_carbon_chemistry(cc, DIC, Alk, T, S) = cc(; DIC, Alk, T, S)
@inline call_carbon_chemistry(cc, DIC, Alk, T, S, silicon, phosphate) = cc(; DIC, Alk, T, S, silicon, phosphate)

# default values from Dickson et al., 2007
@kwdef struct PolynomialVirialCoefficientForCarbonDioxide{FT}
    a :: FT = -1636.75
    b :: FT =  12.0408
    c :: FT = -3.27957 * 10^-2
    d :: FT =  3.16528 * 10^-5
end

@inline (fvc::PolynomialVirialCoefficientForCarbonDioxide)(Tk) = 
    (fvc.a + fvc.b * Tk + fvc.c * Tk^2 + fvc.d * Tk^3) * 10^-6 # cm² / mol to m² / mol

# default values from Dickson et al., 2007
@kwdef struct CrossVirialCoefficientForCarbonDioxide{FT}
    a :: FT =  57.7
    b :: FT =  -0.118
end

@inline (cvc::CrossVirialCoefficientForCarbonDioxide)(Tk) = (cvc.a + cvc.b * Tk) * 10^-6  # cm² / mol to m² / mol

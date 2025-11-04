# default values from Dickson et al., 2007
@kwdef struct PolynomialVirialCoefficientForCarbonDioxide{FT}
    a :: FT = -1636.75
    b :: FT =  12.0408
    c :: FT = -3.27957 * 10^-2
    d :: FT =  3.16528 * 10^-5
end

@inline (fvc::PolynomialVirialCoefficientForCarbonDioxide)(Tk::FT) where FT = 
    (fvc.a + fvc.b * Tk + fvc.c * Tk^2 + fvc.d * Tk^3) * convert(FT, 10^-6) # cm² / mol to m² / mol

# default values from Dickson et al., 2007
@kwdef struct CrossVirialCoefficientForCarbonDioxide{FT}
    a :: FT =  57.7
    b :: FT =  -0.118
end

@inline (cvc::CrossVirialCoefficientForCarbonDioxide)(Tk::FT) where FT = (cvc.a + cvc.b * Tk) * convert(FT, 10^-6)  # cm² / mol to m² / mol
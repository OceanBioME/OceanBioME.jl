"""
    PolynomialSchmidtNumber

Returns schmidt number at a given temperature parameterised in the form
``a + bT + cT² + dT³ + eT⁴``

Form from Wanninkhof, 2014.
"""
struct PolynomialSchmidtNumber{FT}
    a :: FT
    b :: FT
    c :: FT
    d :: FT
    e :: FT
end

@inline (sc::PolynomialSchmidtNumber)(T) = 
    sc.a + sc.b * T + sc.c * T^2 + sc.d * T^3 + sc.e * T^4

"""
    CarbonDioxidePolynomialSchmidtNumber(; a = 2116.8, b = -136.25, c = 4.7353, d = -0.092307, e = 0.0007555)

Schmidt number parameterisation Wanninkhof, 2014 for sea water
"""    
CarbonDioxidePolynomialSchmidtNumber(; a = 2116.8, b = -136.25, c = 4.7353, d = -0.092307, e = 0.0007555) =
    PolynomialSchmidtNumber(a, b, c, d, e)

"""
    OxygenPolynomialSchmidtNumber(; a = 1953.4, b = - 128.0, c = 3.9918, d = -0.050091)

Schmidt number parameterisation Wanninkhof, 2014 for sea water
"""
OxygenPolynomialSchmidtNumber(; a = 1920.4, b = -135.6, c = 5.2122, d = -0.10939, e = 0.00093777) =
    PolynomialSchmidtNumber(a, b, c, d, e)
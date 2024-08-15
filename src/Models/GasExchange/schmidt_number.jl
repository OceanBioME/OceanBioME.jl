
"""
    CarbonDioxidePolynomialSchmidtNumber(; a = 2116.8, b = -136.25, c = 4.7353, d = -0.092307, e = 0.0007555)

Schmidt number parameterisation Wanninkhof, 2014 for sea water
"""    
CarbonDioxidePolynomialSchmidtNumber(; a = 2116.8, b = -136.25, c = 4.7353, d = -0.092307, e = 0.0007555) =
    PolynomialParameterisation{4}((a, b, c, d, e))

"""
    OxygenPolynomialSchmidtNumber(; a = 1953.4, b = - 128.0, c = 3.9918, d = -0.050091)

Schmidt number parameterisation Wanninkhof, 2014 for sea water
"""
OxygenPolynomialSchmidtNumber(; a = 1920.4, b = -135.6, c = 5.2122, d = -0.10939, e = 0.00093777) =
    PolynomialParameterisation{4}((a, b, c, d, e))
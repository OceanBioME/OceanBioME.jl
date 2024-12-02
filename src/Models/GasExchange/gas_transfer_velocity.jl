module ScaledGasTransferVelocity

export SchmidtScaledTransferVelocity

using Oceananigans.Units, Adapt

using OceanBioME.Models.GasExchangeModel: PolynomialParameterisation

import Adapt: adapt_structure

"""
    SchmidtScaledTransferVelocity(; schmidt_number, 
                                    base_transfer_velocity = Ho06())

Returns a model for gas transfer velocity which depends on the `u₁₀`, the 10m-wind, and 
`T`emperature. The model is of the typical form:

    k(u₁₀, T) = k₆₆₀(u₁₀) √(660/Sc(T))

The `base_transfer_velocity` (k₆₆₀) is typically an empirically derived gas transfer velocity
normalised by the Scmidt number for CO₂ at 20°C (660), and the `schmidt_number` (Sc) is a parameterisation
of the gas specific Schmidt number.
"""
@kwdef struct SchmidtScaledTransferVelocity{KB, SC} 
  base_transfer_velocity :: KB = Ho06()
          schmidt_number :: SC
end

(k::SchmidtScaledTransferVelocity)(u₁₀, T) = k.base_transfer_velocity(u₁₀) * (k.schmidt_number(T) / 660)^(-1/2)

Adapt.adapt_structure(to, k::SchmidtScaledTransferVelocity) = SchmidtScaledTransferVelocity(adapt(to, k.base_transfer_velocity),
                                                                                            adapt(to, k.schmidt_number))

summary(::SchmidtScaledTransferVelocity{KB, SC}) where {KB, SC} = "SchmidtScaledTransferVelocity{$(nameof(KB)), $(nameof(SC))}"
show(io::IO, k::SchmidtScaledTransferVelocity{KB, SC}) where {KB, SC} = 
    println(io, summary(k), "\n",
                "    k = k₆₆₀(u₁₀) √(660/Sc(T)),\n",
                "    Sc(T): $(nameof(SC)),\n",
                "    k₆₆₀(u₁₀) : $(nameof(KB))")

"""
    Wanninkhof99(FT = Float64; scale_factor = 0.0283 / hour / 100)

Cubic k₆₆₀ parameterisation of Wanninkhof & McGillis (1999) suitable for 
short term, in situ wind products.
"""
Wanninkhof99(FT = Float64; scale_factor = 0.0283 / hour / 100) = PolynomialParameterisation{3}(FT; coefficients = (0, 0, 0, scale_factor))

"""
    Ho06(FT = Float64; scale_factor = 0.266 / hour / 100)

Quadratic k₆₆₀ parameterisation of Ho et al. (2006) suitable for the QuickSCAT satellite and short-term 
steady wind product.
"""
Ho06(FT = Float64; scale_factor = 0.266 / hour / 100) = PolynomialParameterisation{2}(FT; coefficients = (0, 0, scale_factor))

"""
   Nightingale00(FT = Float64; linear = 0.333 / hour / 100, quadratic = 0.222 / hour / 100)

Cubic k₆₆₀ parameterisation of Nightingale et al. (2000) suitable for 
short term, in situ wind products (?).
"""
Nightingale00(FT = Float64; linear = 0.333 / hour / 100, quadratic = 0.222 / hour / 100) =
    PolynomialParameterisation{2}(FT; coefficients = (0, linear, quadratic))

"""
    McGillis01(FT = Float64; constant = 3.3 / hour / 100, cubic = 0.026 / hour / 100)

Cubic k₆₆₀ parameterisation of McGillis et al. (2001) suitable for 
short term, in situ wind products.
"""
McGillis01(FT = Float64; constant = 3.3 / hour / 100, cubic = 0.026 / hour / 100) =
    PolynomialParameterisation{3}(FT; coefficients = (constant, 0, 0, cubic))

"""
    Sweeny07(FT = Float64; scale_factor = 0.27 / hour / 100)

Quadratic k₆₆₀ parameterisation of Sweeny et al. (2007) suitable for the
NCEP/NCAR reanalysis 1 product
"""
Sweeny07(FT = Float64; scale_factor = 0.27 / hour / 100) = PolynomialParameterisation{2}(FT; coefficients = (0, 0, scale_factor))

"""
    Wanninkhof09(FT = Float64; constant = 3 / hour / 100, linear = 0.1 / hour / 100, quadratic = 0.064 / hour / 100, cubic = 0.011 / hour / 100)

Cubic k₆₆₀ parameterisation of Wanninkhof et al (2009) suitable for the
Cross-Calibrated Multi-Platform (CCMP) Winds product
"""
Wanninkhof09(FT = Float64; constant = 3 / hour / 100, linear = 0.1 / hour / 100, quadratic = 0.064 / hour / 100, cubic = 0.011 / hour / 100) =
    PolynomialParameterisation{3}(FT; coefficients = (constant, linear, quadratic, cubic))


"""
    Wanninkhof14(FT = Float64; scale_factor = 0.251 / hour / 100)

Quadratic k₆₆₀ parameterisation of Wanninkhof et al (2014) suitable for the
Cross-Calibrated Multi-Platform (CCMP) Winds product
"""
Wanninkhof14(FT = Float64; scale_factor = 0.251 / hour / 100) = PolynomialParameterisation{2}(FT; coefficients = (0, 0, scale_factor))

"""
    ERA5(FT = Float64; scale_factor = 0.270875 / hour / 100)

Quadratic k₆₆₀ parameterisation calibrated to give 16.5 cm/hr global average (reccomended Naegler, 2009)
for the ERA5 wind product by SeaFlux/Luke Gregor et al. (2023).
"""
ERA5(FT = Float64; scale_factor = 0.270875 / hour / 100) = PolynomialParameterisation{2}(FT; coefficients = (0, 0, scale_factor))

"""
    JRA55(FT = Float64; scale_factor = 0.2601975 / hour / 100)

Quadratic k₆₆₀ parameterisation calibrated to give 16.5 cm/hr global average (reccomended Naegler, 2009)
for the JRA55 wind product by SeaFlux/Luke Gregor et al. (2023).
"""
JRA55(FT = Float64; scale_factor = 0.2601975 / hour / 100) = PolynomialParameterisation{2}(FT; coefficients = (0, 0, scale_factor))

"""
    NCEP1(FT = Float64; scale_factor = 0.2866424 / hour / 100)

Quadratic k₆₆₀ parameterisation calibrated to give 16.5 cm/hr global average (reccomended Naegler, 2009)
for the NCEP1 wind product by SeaFlux/Luke Gregor et al. (2023).
"""
NCEP1(FT = Float64; scale_factor = 0.2866424 / hour / 100) = PolynomialParameterisation{2}(FT; coefficients = (0, 0, scale_factor))

"""
    CCMP2(FT = Float64; scale_factor = 0.256789 / hour / 100)

Quadratic k₆₆₀ parameterisation calibrated to give 16.5 cm/hr global average (reccomended Naegler, 2009)
for the CCMP2 wind product by SeaFlux/Luke Gregor et al. (2023).
"""
CCMP2(FT = Float64; scale_factor = 0.256789 / hour / 100) = PolynomialParameterisation{2}(FT; coefficients = (0, 0, scale_factor))

end # module
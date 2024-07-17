
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

"""
    K0(; constant = -162.8301
         inverse_T =  218.2968 * 100
         log_T =  90.9241
         T² = -1.47696 / 100^2
         S =  0.025695
         ST = -0.025225 / 100
         ST² =  0.0049867 / 100^2)

Parameterisation for carbon dioxide solubility equilibrium constant.

    CO₂(g) ⇌ CO₂*(aq)

    K₀ = [CO₂*(aq)]/f(CO₂)

Default values from Weiss & Price (1980, Mar. Chem., 8, 347-359; Eq 13 with table 6 values).
"""
@kwdef struct K0{FT}
     constant :: FT = -162.8301
    inverse_T :: FT =  218.2968 * 100
        log_T :: FT =  90.9241
           T² :: FT = -1.47696 / 100^2
            S :: FT =  0.025695
           ST :: FT = -0.025225 / 100
          ST² :: FT =  0.0049867 / 100^2
end

@inline (c::K0)(T, S) = exp(c.constant 
                            + c.inverse_T / T
                            + c.log_T * (log(T) - log(100))
                            + c.T² * T^2
                            + (c.S + c.ST * T + c.ST² * T^2) * S)

"""
    K1(; constant =  62.008, 
         inverse_T = -3670.7,
         log_T = -9.7944,
         S =  0.0118,
         S² = -0.000116)

Parameterisation for aquious carbon dioxide - bicarbonate dissociation equilibrium constant.

    CO₂*(aq) + H₂O ⇌ H₂CO₃ ⇌ HCO₃⁻ + H⁺

    K₁ = [H⁺][HCO₃⁻]/[CO₂*]

Default values from Lueke, et. al (2000, Mar. Chem., 70, 105–119;).
"""
@kwdef struct K1{FT}
     constant :: FT =  62.008
    inverse_T :: FT = -3670.7
        log_T :: FT = -9.7944
            S :: FT =  0.0118
           S² :: FT = -0.000116
end

@inline (c::K1)(T, S) = 10 ^ (c.constant + c.inverse_T / T + c.log_T * log(T) + c.S * S + c.S² * S^2)

"""
    K2(; constant =  62.008, 
         inverse_T = -3670.7,
         log_T = -9.7944,
         S =  0.0118,
         S² = -0.000116)

Parameterisation for bicarbonate dissociation equilibrium constant.

    HCO₃⁻ ⇌ CO₃²⁻ + H⁺

    K₂ = [H⁺][CO₃²⁻]/[HCO₃⁻]

Default values from Lueke, et. al (2000, Mar. Chem., 70, 105–119).
"""
@kwdef struct K2{FT}
     constant :: FT = -4.777
    inverse_T :: FT = -1394.7
            S :: FT =  0.0184
           S² :: FT = -0.000118
        log_T :: FT =  0.0
end

@inline (c::K2)(T, S) = 10 ^ (c.constant + c.inverse_T / T + c.S * S + c.S² * S^2 + c.log_T * log(T))
    
"""
    KB(; constant =  148.0248,
         inverse_T = -8966.90,
         invsese_T_sqrt_S = -2890.53,
         invsese_T_S = -77.942,
         invsese_T_sqrt_S³ =  1.728,
         inverse_T_S² = -0.0996,
         sqrt_S = 137.1942,
         S = 1.62142,
         log_T = -24.4344,
         log_T_sqrt_S = -25.085,
         S_log_T = -0.2474,
         T_sqrt_S =  0.053105)

Parameterisation for boric acid dissociation equilibrium constant.

    B(OH)₃ + H₂O ⇌ B(OH)₄⁻ + H⁺

    Kᵇ = [H⁺][B(OH)₄⁻]/[B(OH)₃]

Default values from Dickson (1990, Deep-Sea Res., 37, 755–766).
"""
@kwdef struct KB{FT}
          constant :: FT =  148.0248
         inverse_T :: FT = -8966.90
  invsese_T_sqrt_S :: FT = -2890.53
       invsese_T_S :: FT = -77.942
 invsese_T_sqrt_S³ :: FT =  1.728
      inverse_T_S² :: FT = -0.0996
            sqrt_S :: FT = 137.1942
                 S :: FT = 1.62142
             log_T :: FT = -24.4344
      log_T_sqrt_S :: FT = -25.085
           S_log_T :: FT = -0.2474
          T_sqrt_S :: FT =  0.053105
end

@inline (c::KB)(T, S) = exp(c.constant 
                            + (c.inverse_T + c.invsese_T_sqrt_S * √S + c.invsese_T_S * S + c.invsese_T_sqrt_S³ * S^1.5 + c.inverse_T_S² * S^2) / T
                            + c.sqrt_S * √S
                            + c.S * S
                            + (c.log_T + c.log_T_sqrt_S * √S + c.S_log_T * S ) * log(T)
                            + c.T_sqrt_S * √S * T)
 
"""
    KW(; constant =  148.9652,
         inverse_T = -13847.26,
         log_T = -23.6521,
         sqrt_S = -5.977,
         inverse_T_sqrt_S =  118.67,
         log_T_sqrt_S =  1.0495,
         S = -0.01615)

Parameterisation for water dissociation equilibrium constant.

    H₂O ⇌ OH⁻ + H⁺

    Kʷ = [H⁺][OH⁻]

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
@kwdef struct KW{FT}
          constant :: FT =  148.9652
         inverse_T :: FT = -13847.26
             log_T :: FT = -23.6521
            sqrt_S :: FT = -5.977
  inverse_T_sqrt_S :: FT =  118.67
      log_T_sqrt_S :: FT =  1.0495
                 S :: FT = -0.01615
end

@inline (c::KW)(T, S) = exp(c.constant
                            + c.inverse_T / T
                            + c.log_T * log(T)
                            + (c.sqrt_S + c.inverse_T_sqrt_S / T + c.log_T_sqrt_S * log(T))* √S
                            + c.S * S)

"""
    IonicStrength(; a =  19.924,
                    b =  1000.0,
                    c = -1.005)

Parameterisation of the ionic strength of sea water.

    Is(S) = aS / (b + cS)

Default values from Dickson (1990, Chem. Thermodyn., 22, 113–127).
"""
@kwdef struct IonicStrength{FT}
    a :: FT =  19.924
    b :: FT =  1000.0
    c :: FT = -1.005
end

@inline (c::IonicStrength)(S) = c.a * S / (c.b + c.c * S)

"""
    KS(; constant =  148.9652,
         inverse_T = -13847.26,
         log_T = -23.6521,
         sqrt_S = -5.977,
         inverse_T_sqrt_S =  118.67,
         log_T_sqrt_S =  1.0495,
         S = -0.01615)

Parameterisation for bisulfate dissociation equilibrium constant.

    HSO₄⁻ ⇌ SO₄⁻ + H⁺

    Kˢ = [H⁺][SO₄⁻]/[HSO₄⁻] 

Default values from Dickson (1990, Chem. Thermodyn., 22, 113–127).
"""
@kwdef struct KS{IS, FT}
    ionic_strength :: IS = IonicStrength()

          constant :: FT = 141.328
         inverse_T :: FT = -4276.1
             log_T :: FT = -23.093
            sqrt_S :: FT =  324.57
  inverse_T_sqrt_S :: FT = -13856.0
      log_T_sqrt_S :: FT = -47.986
                Is :: FT = -771.54
      inverse_T_Is ::FT =  35474.0
          log_T_Is :: FT =  114.723
 inverse_T_sqrt_S³ :: FT =  2698.0
      inverse_T_S² :: FT =  1776.0
             log_S :: FT = -0.001005
end

@inline (c::KS)(T, S, Is = c.ionic_strength(S)) = exp(c.constant
                                                      + c.inverse_T / T
                                                      + c.log_T * log(T)
                                                      + (c.sqrt_S + c.inverse_T_sqrt_S / T + c.log_T_sqrt_S * log(T)) * √S
                                                      + (c.Is + c.inverse_T_Is / T + c.log_T_Is * log(T)) * Is
                                                      + c.inverse_T_sqrt_S³ * S^1.5 / T
                                                      + c.inverse_T_S² * S^2 / T
                                                      + log(1 + c.log_S * S))

"""
    KF(; ionic_strength = IonicStrength(),
          sulfate_constant = KS(; ionic_strength),
          constant = -12.641,
          inverse_T =  1590.2,
          sqrt_S =  1.525,
          log_S = -0.001005,
          log_S_KS = 0.1400 / 96.062 / 1.80655)

Parameterisation for hydrogen fluoride dissociation equilibrium constant.

    HF ⇌ F⁻ + H⁺

    Kᶠ = [H⁺][F⁻]/[HF] 

Default values from Perez and Fraga (1987, Mar. Chem., 21, 161–168).
"""
@kwdef struct KF{IS, KS, FT}
     ionic_strength :: IS = IonicStrength()
   sulfate_constant :: KS = KS(; ionic_strength)

           constant :: FT = -12.641
          inverse_T :: FT =  1590.2
             sqrt_S :: FT =  1.525
              log_S :: FT = -0.001005
           log_S_KS :: FT = 0.1400 / 96.062 / 1.80655
end

@inline (c::KF)(T, S, Is = c.ionic_strength(S), KS = c.sulfate_constant(T, S, Is)) =
    exp(c.constant
        + c.inverse_T / T
        + c.sqrt_S * √S
        + log(1 + c.log_S * S)
        + log(1 + c.log_S_KS * S / KS))

"""
    KP(constant,
       inverse_T,
       log_T,
       sqrt_S,
       inverse_T_sqrt_S,
       S,
       inverse_T_S)

Generic equilibrium constant parameterisation of the form used by 
Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677) for phosphoric 
acid dissociation.
"""
struct KP{FT}
          constant :: FT
         inverse_T :: FT
             log_T :: FT
            sqrt_S :: FT
  inverse_T_sqrt_S :: FT
                 S :: FT
       inverse_T_S :: FT
end

@inline (c::KP)(T, S) = exp(c.constant
                            + c.inverse_T / T
                            + c.log_T * log(T)
                            + (c.sqrt_S + c.inverse_T_sqrt_S / T) * √S
                            + (c.S + c.inverse_T_S / T) * S)

"""
    KP1(; constant = 115.525,
          inverse_T = -4576.752,
          log_T = - 18.453,
          sqrt_S = 0.69171,
          inverse_T_sqrt_S = -106.736,
          S = -0.01844,
          inverse_T_S = -0.65643)

Instance of `KP` returning the first phosphocic acid equilibrium constant.

    H₃PO₄ ⇌ H₂PO₄⁻ + H⁺

    Kᵖ¹ = [H⁺][H₂PO₄]/[H₃PO₄] 

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
KP1(; constant = 115.525,
      inverse_T = -4576.752,
      log_T = - 18.453,
      sqrt_S = 0.69171,
      inverse_T_sqrt_S = -106.736,
      S = -0.01844,
      inverse_T_S = -0.65643) = KP(constant,
                                   inverse_T,
                                   log_T,
                                   sqrt_S,
                                   inverse_T_sqrt_S,
                                   S,
                                   inverse_T_S)

"""
    KP2(; constant = 172.0883,
          inverse_T = -8814.715,
          log_T = -27.927,
          sqrt_S = 1.3566,
          inverse_T_sqrt_S = -160.340,
          S = -0.05778,
          inverse_T_S = 0.37335)

Instance of `KP` returning the second phosphocic acid equilibrium constant.

    H₂PO₄⁻ ⇌ HPO₄ + H⁺

    Kᵖ² = [H⁺][HPO₄]/[H₂PO₄⁻] 

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
KP2(; constant = 172.0883,
      inverse_T = -8814.715,
      log_T = -27.927,
      sqrt_S = 1.3566,
      inverse_T_sqrt_S = -160.340,
      S = -0.05778,
      inverse_T_S = 0.37335) = KP(constant,
                                   inverse_T,
                                   log_T,
                                   sqrt_S,
                                   inverse_T_sqrt_S,
                                   S,
                                   inverse_T_S)

"""
    KP3(; constant = - 18.141,
          inverse_T = -3070.75,
          log_T = 0.0,
          sqrt_S = 2.81197,
          inverse_T_sqrt_S = 17.27039,
          S = -0.09984,
          inverse_T_S = -44.99486)

Instance of `KP` returning the third phosphocic acid equilibrium constant.

    HPO₄⁻ ⇌ PO₄ + H⁺

    Kᵖ³ = [H⁺][PO₄]/[HPO₄⁻] 

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
KP3(; constant = - 18.141,
      inverse_T = -3070.75,
      log_T = 0.0,
      sqrt_S = 2.81197,
      inverse_T_sqrt_S = 17.27039,
      S = -0.09984,
      inverse_T_S = -44.99486) = KP(constant,
                                   inverse_T,
                                   log_T,
                                   sqrt_S,
                                   inverse_T_sqrt_S,
                                   S,
                                   inverse_T_S)

"""
    KSi(; ionic_strength = IonicStrength(),
          constant =  117.385,
          inverse_T = -8904.2,
          log_T = -19.334,
          sqrt_Is = 3.5913,
          inverse_T_sqrt_Is = -458.79,
          Is = -1.5998,
          inverse_T_Is = 188.74,
          Is² = 0.07871,
          inverse_T_Is² = -12.1652,
          log_S = -0.001005)

Parameterisation for silicic acid dissociation equilibrium constant.

    Si(OH)₄ ⇌ SiO(OH)₃⁻ + H⁺

    Kʷ = [H⁺][SiO(OH)₃⁻]/[Si(OH)₄]

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
@kwdef struct KSi{IS, FT}
    ionic_strength :: IS = IonicStrength()

    constant :: FT =  117.385
    inverse_T :: FT = -8904.2
    log_T :: FT = -19.334
    sqrt_Is :: FT = 3.5913
    inverse_T_sqrt_Is :: FT = -458.79
    Is :: FT = -1.5998
    inverse_T_Is :: FT = 188.74
    Is² :: FT = 0.07871
    inverse_T_Is² :: FT = -12.1652
    log_S :: FT = -0.001005
end

@inline (c::KSi)(T, S, Is = c.ionic_strength(S)) = exp(c.constant
                                                       + c.inverse_T / T
                                                       + c.log_T * log(T)
                                                       + (c.sqrt_Is + c.inverse_T_sqrt_Is / T) * √Is
                                                       + (c.Is + c.inverse_T_Is / T) * Is
                                                       + (c.Is² + c.inverse_T_Is² / T) * Is^2
                                                       + log(1 + c.log_S * S))

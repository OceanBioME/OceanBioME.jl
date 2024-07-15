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

@kwdef struct K1{FT}
     constant :: FT =  62.008
    inverse_T :: FT = -3670.7
        log_T :: FT = -9.7944
            S :: FT =  0.0118
           S² :: FT = -0.000116
end

@inline (c::K1)(T, S) = 10 ^ (c.constant + c.inverse_T / T + c.log_T * log(T) + c.S * S + c.S² * S^2)

@kwdef struct K2{FT}
     constant :: FT = -4.777
    inverse_T :: FT = -1394.7
            S :: FT =  0.0184
           S² :: FT = -0.000118
        log_T :: FT =  0.0
end

@inline (c::K2)(T, S) = 10 ^ (c.constant + c.inverse_T / T + c.S * S + c.S² * S^2 + c.log_T * log(T))
    
@kwdef struct KB{FT}
          constant :: FT =  148.0248
         inverse_T :: FT = -8966.90
  invsese_T_sqrt_S :: FT = - 2890.53
       invsese_T_S :: FT = - 77.942
 invsese_T_sqrt_S³ :: FT =  1.728
      inverse_T_S² :: FT = - 0.0996
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

@kwdef struct KW{FT}
          constant :: FT =  148.9652
         inverse_T :: FT = -13847.26
             log_T :: FT = -23.6521
            sqrt_S :: FT = -5.977
  inverse_T_sqrt_S :: FT =  118.67
      log_T_sqrt_S :: FT =  1.0495
                 S :: FT = - 0.01615
end

@inline (c::KW)(T, S) = exp(c.constant
                            + c.inverse_T / T
                            + c.log_T * log(T)
                            + (c.sqrt_S + c.inverse_T_sqrt_S / T + c.log_T_sqrt_S * log(T))* √S
                            + c.S * S)

@kwdef struct IonicStrength{FT}
    a :: FT =  19.924
    b :: FT =  1000.0
    c :: FT = -1.005
end

@inline (c::IonicStrength)(S) = c.a * S / (c.b + c.c * S)

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

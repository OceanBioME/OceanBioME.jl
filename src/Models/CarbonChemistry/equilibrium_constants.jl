"""
    PressureCorrection(FT=Float64;
                       a₀, a₁, a₂,
                       b₀, b₁, b₂,
                       R = 83.14472)

Parameterisation for the pressure effect on thermodynamic constants.

Form from Millero, F. J. (2007, Chemical Reviews, 107(2), 308–341).
"""
struct PressureCorrection{FT}
    a₀ :: FT
    a₁ :: FT
    a₂ :: FT
    b₀ :: FT
    b₁ :: FT

    R  :: FT 

    function PressureCorrection(FT = Float64; 
                                a₀, a₁, a₂, 
                                b₀, b₁,
                                R = 83.14472)

        return new{FT}(a₀, a₁, a₂, b₀, b₁, R)
    end
end

@inline function (pc::PressureCorrection)(Tk, P)
    Tc = Tk - 273.15

    ΔV = pc.a₀ + pc.a₁ * Tc + pc.a₂ * Tc^2
    Δκ = pc.b₀ + pc.b₁ * Tc

    RT = pc.R * Tk

    return exp((-ΔV + 0.5 * Δκ * P) * P / RT)
end

@inline (pc::PressureCorrection{FT})(T, ::Nothing) where FT = convert(FT, 1)

summary(::IO, ::PressureCorrection) = string("Equilibrium constant pressure correction")
show(io::IO, pc::PressureCorrection) = print(io, "Equilibrium constant pressure correction\n",
    "    ΔV = $(pc.a₀) + $(pc.a₁) Tc + $(pc.a₂) * Tc²\n",
    "    Δκ = $(pc.b₀) + $(pc.b₁) T\n",
    "    ln(kᵖ/k⁰) = (-ΔV + Δκ P / 2) * P / RT")

"""
    K0(; constant = -60.2409,
         inverse_T =  93.4517 * 100,
         log_T =  23.3585,
         T² =  0.0,
         S =  0.023517,
         ST = -0.023656 / 100,
         ST² =  0.0047036 / 100^2)

Parameterisation for carbon dioxide solubility equilibrium constant.

    CO₂(g) ⇌ CO₂*(aq)

    K₀ = [CO₂*(aq)]/f(CO₂)

Default values from Weiss, R.F. (1974, Mar. Chem., 2, 203–215).
"""
@kwdef struct K0{FT}
     constant :: FT = -60.2409
    inverse_T :: FT =  93.4517 * 100
        log_T :: FT =  23.3585
           T² :: FT =  0.0
            S :: FT =  0.023517
           ST :: FT = -0.023656 / 100
          ST² :: FT =  0.0047036 / 100^2
end

@inline (c::K0)(T::FT, S; P = nothing) where FT = 
    exp(c.constant 
        + c.inverse_T / T
        + c.log_T * (log(T) - log(convert(FT, 100)))
        + c.T² * T^convert(FT, 2)
        + (c.S + c.ST * T + c.ST² * T^convert(FT, 2)) * S)

summary(::IO, ::K0) = string("Solubility constant")
show(io::IO, k0::K0) = print(io, "Solubility constant\n",
    "    ln(k₀/k°) = $(k0.constant) + $(k0.inverse_T) / T + $(k0.log_T) (log(T) - log(100)) + $(k0.T²) T² + ($(k0.S) + $(k0.ST) T + $(k0.ST²) T²)S")

"""
    K1(FT = Float64;
       constant =  61.2172,
       inverse_T = -3633.86,
       log_T = -9.67770,
       S =  0.011555,
       S² = -0.0001152,
       pressure_correction = PressureCorrection(FT; a₀=-25.50, a₁=0.1271, a₂=0.0, b₀=-0.00308, b₁=0.0000877))

Parameterisation for aquious carbon dioxide - bicarbonate dissociation equilibrium constant.

    CO₂*(aq) + H₂O ⇌ H₂CO₃ ⇌ HCO₃⁻ + H⁺

    K₁ = [H⁺][HCO₃⁻]/[CO₂*]

Default values from Lueker et al. (2000, Mar. Chem., 70: 105–119).
"""
struct K1{FT, PC}
               constant :: FT
              inverse_T :: FT
                  log_T :: FT
                      S :: FT
                     S² :: FT

    pressure_correction :: PC

    function K1(FT = Float64;
                constant =  61.2172,
                inverse_T = -3633.86,
                log_T = -9.67770,
                S =  0.011555,
                S² = -0.0001152,
                pressure_correction::PC = PressureCorrection(FT; a₀=-25.50, a₁=0.1271, a₂=0.0, b₀=-0.00308, b₁=0.0000877)) where PC
        
        return new{FT, PC}(constant, inverse_T, log_T, S, S², pressure_correction)
    end
end

@inline (c::K1)(T, S; P = nothing) = 
    c.pressure_correction(T, P) *
    10 ^ (c.constant + c.inverse_T / T + c.log_T * log(T) + c.S * S + c.S² * S^2)

summary(::IO, ::K1) = string("First carbon dioxide dissociation constant")
show(io::IO, k1::K1) = print(io, "First carbon dioxide dissociation constant\n",
    "    log₁₀(k₁/k°) = $(k1.constant) + $(k1.inverse_T) / T + $(k1.log_T) log(T) + $(k1.S) S + $(k1.S²) S²")

"""
    K2(FT = Float64;
       constant = -25.9290,
       inverse_T = -471.78,
       log_T = 3.16967,
       S = 0.01781,
       S² = -0.0001122,
       pressure_correction = PressureCorrection(FT; a₀=-15.82, a₁=-0.0219, a₂=0.0, b₀=0.00113, b₁=-0.0001475))

Parameterisation for bicarbonate dissociation equilibrium constant.

    HCO₃⁻ ⇌ CO₃²⁻ + H⁺

    K₂ = [H⁺][CO₃²⁻]/[HCO₃⁻]

Default values from Lueker et al. (2000, Mar. Chem., 70: 105–119).
"""
struct K2{FT, PC}
               constant :: FT
              inverse_T :: FT
                  log_T :: FT
                      S :: FT
                     S² :: FT

    pressure_correction :: PC

    function K2(FT = Float64;
                constant = -25.9290,
                inverse_T = -471.78,
                log_T = 3.16967,
                S = 0.01781,
                S² = -0.0001122,
                pressure_correction::PC = PressureCorrection(FT; a₀=-15.82, a₁=-0.0219, a₂=0.0, b₀=0.00113, b₁=-0.0001475)) where PC
        
        return new{FT, PC}(constant, inverse_T, log_T, S, S², pressure_correction)
    end
end

@inline (c::K2)(T, S; P = nothing) = 
    c.pressure_correction(T, P) *
    10 ^ (c.constant + c.inverse_T / T + c.S * S + c.S² * S^2 + c.log_T * log(T))

summary(::IO, ::K2) = string("Second carbon dioxide dissociation constant")
show(io::IO, k2::K2) = print(io, "Second carbon dioxide dissociation constant\n",
    "    log₁₀(k₂/k°) = $(k2.constant) + $(k2.inverse_T) / T + $(k2.log_T) log(T) + $(k2.S) S + $(k2.S²) S²")

"""
    KB(FT = Float64;
       constant =  148.0248,
       inverse_T = -8966.90,
       inverse_T_sqrt_S = -2890.53,
       inverse_T_S = -77.942,
       inverse_T_sqrt_S³ =  1.728,
       inverse_T_S² = -0.0996,
       sqrt_S = 137.1942,
       S = 1.62142,
       log_T = -24.4344,
       log_T_sqrt_S = -25.085,
       S_log_T = -0.2474,
       T_sqrt_S =  0.053105,
       pressure_correction = PressureCorrection(FT; a₀=-29.48, a₁=0.1622, a₂=-0.0026080, b₀=-0.00284, b₁=0.0))

Parameterisation for boric acid equilibrium with water.

    B(OH)₃ + H₂O ⇌ B(OH)₄⁻ + H⁺

    Kᵇ = [H⁺][B(OH)₄⁻]/[B(OH)₃]

Default values from Dickson (1990, Deep-Sea Res., 37, 755–766).
"""
struct KB{FT, PC}
               constant :: FT
              inverse_T :: FT
       inverse_T_sqrt_S :: FT
            inverse_T_S :: FT
      inverse_T_sqrt_S³ :: FT
           inverse_T_S² :: FT
                 sqrt_S :: FT
                      S :: FT
                  log_T :: FT
           log_T_sqrt_S :: FT
                S_log_T :: FT
               T_sqrt_S :: FT

    pressure_correction :: PC

    function KB(FT = Float64;
                constant =  148.0248,
                inverse_T = -8966.90,
                inverse_T_sqrt_S = -2890.53,
                inverse_T_S = -77.942,
                inverse_T_sqrt_S³ =  1.728,
                inverse_T_S² = -0.0996,
                sqrt_S = 137.1942,
                S = 1.62142,
                log_T = -24.4344,
                log_T_sqrt_S = -25.085,
                S_log_T = -0.2474,
                T_sqrt_S =  0.053105,
                pressure_correction::PC = PressureCorrection(FT; a₀=-29.48, a₁=0.1622, a₂=-0.0026080, b₀=-0.00284, b₁=0.0)) where PC
        
        return new{FT, PC}(constant, inverse_T, 
                           inverse_T_sqrt_S, inverse_T_S, 
                           inverse_T_sqrt_S³, inverse_T_S², 
                           sqrt_S, S, 
                           log_T, log_T_sqrt_S, 
                           S_log_T, T_sqrt_S, 
                           pressure_correction)
    end
end

@inline (c::KB)(T::FT, S; P = nothing) where FT = 
    c.pressure_correction(T, P) *
    exp(c.constant 
        + (c.inverse_T + c.inverse_T_sqrt_S * √S + c.inverse_T_S * S + c.inverse_T_sqrt_S³ * S^convert(FT, 1.5) + c.inverse_T_S² * S^convert(FT, 2)) / T
        + c.sqrt_S * √S
        + c.S * S
        + (c.log_T + c.log_T_sqrt_S * √S + c.S_log_T * S ) * log(T)
        + c.T_sqrt_S * √S * T)

summary(::IO, ::KB) = string("Boric acid dissociation constant")
show(io::IO, c::KB) = print(io, "Boric acid dissociation constant\n",
    "    ln(kᵇ/k°) = $(c.constant) + ($(c.inverse_T) + $(c.inverse_T_sqrt_S) √S + $(c.inverse_T_S) S + $(c.inverse_T_sqrt_S³) √S³ + $(c.inverse_T_S²) S²) / T
                + $(c.sqrt_S) * √S
                + $(c.S) * S
                + ($(c.log_T) + $(c.log_T_sqrt_S) √S + $(c.S_log_T) S ) * log(T)
                + $(c.T_sqrt_S) * √S * T")

"""
    KW(; constant =  148.9652,
         inverse_T = -13847.26,
         log_T = -23.6521,
         sqrt_S = -5.977,
         inverse_T_sqrt_S =  118.67,
         log_T_sqrt_S =  1.0495,
         S = -0.01615,
         pressure_correction = 
            PressureCorrection(; a₀=-20.02, a₁=0.1119, a₂=-0.001409, b₀=-0.00513, b₁=0.0000794))

Parameterisation for water dissociation equilibrium constant.

    H₂O ⇌ OH⁻ + H⁺

    Kʷ = [H⁺][OH⁻]

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
struct KW{FT, PC}
               constant :: FT 
              inverse_T :: FT 
                  log_T :: FT 
                 sqrt_S :: FT 
       inverse_T_sqrt_S :: FT 
           log_T_sqrt_S :: FT 
                      S :: FT 

    pressure_correction :: PC

    function KW(FT = Float64;
                constant =  148.9652,
                inverse_T = -13847.26,
                log_T = -23.6521,
                sqrt_S = -5.977,
                inverse_T_sqrt_S =  118.67,
                log_T_sqrt_S =  1.0495,
                S = -0.01615,
                pressure_correction::PC = PressureCorrection(FT; a₀=-20.02, a₁=0.1119, a₂=-0.001409, b₀=-0.00513, b₁=0.0000794)) where PC
        
        return new{FT, PC}(constant, inverse_T, 
                           log_T, sqrt_S, 
                           inverse_T_sqrt_S, log_T_sqrt_S, S,
                           pressure_correction)
    end
end

@inline (c::KW)(T, S; P = nothing) =
    c.pressure_correction(T, P) *
    exp(c.constant
        + c.inverse_T / T
        + c.log_T * log(T)
        + (c.sqrt_S + c.inverse_T_sqrt_S / T + c.log_T_sqrt_S * log(T))* √S
        + c.S * S)

summary(::IO, ::KW) = string("Water dissociation constant")
show(io::IO, c::KW) = print(io, "Water dissociation constant\n",
    "    ln(kʷ/k°) = $(c.constant)
                + $(c.inverse_T) / T
                + $(c.log_T) log(T)
                + ($(c.sqrt_S) + $(c.inverse_T_sqrt_S) / T + $(c.log_T_sqrt_S) log(T)) √S
                + $(c.S) * S")


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

summary(::IO, ::IonicStrength) = string("Ionic strength")
show(io::IO, c::IonicStrength) = print(io, "Ionic strength\n",
    "    Is = $(c.a) S/($(c.b) + $(c.c) S)")

"""
    KS(; constant =  148.9652,
         inverse_T = -13847.26,
         log_T = -23.6521,
         sqrt_S = -5.977,
         inverse_T_sqrt_S =  118.67,
         log_T_sqrt_S =  1.0495,
         S = -0.01615,
         pressure_correction = 
            PressureCorrection(; a₀=-18.03, a₁=0.0466, a₂=0.000316, b₀=-0.00453, b₁=0.00009))

Parameterisation for bisulfate dissociation equilibrium constant.

    HSO₄⁻ ⇌ SO₄²⁻ + H⁺

    Kˢ = [H⁺][SO₄²⁻]/[HSO₄⁻] 

Default values from Dickson (1990, Chem. Thermodyn., 22, 113–127).
"""
struct KS{IS, FT, PC}
         ionic_strength :: IS
     
               constant :: FT
              inverse_T :: FT
                  log_T :: FT
                sqrt_Is :: FT
      inverse_T_sqrt_Is :: FT
          log_T_sqrt_Is :: FT
                     Is :: FT
           inverse_T_Is :: FT
               log_T_Is :: FT
     inverse_T_sqrt_Is³ :: FT
          inverse_T_Is² :: FT
                  log_S :: FT

    pressure_correction :: PC

    function KS(FT = Float64;
                ionic_strength::IS = IonicStrength{FT}(),
                constant = 141.328,
                inverse_T = -4276.1,
                log_T = -23.093,
                sqrt_Is =  324.57,
                inverse_T_sqrt_Is = -13856.0,
                log_T_sqrt_Is = -47.986,
                Is = -771.54,
                inverse_T_Is =  35474.0,
                log_T_Is =  114.723,
                inverse_T_sqrt_Is³ = -2698.0,
                inverse_T_Is² =  1776.0,
                log_S = -0.001005,
                pressure_correction::PC =
                    PressureCorrection(FT; a₀=-18.03, a₁=0.0466, a₂=0.000316, b₀=-0.00453, b₁=0.00009)) where {IS, PC}

        return new{IS, FT, PC}(ionic_strength,
                               constant, inverse_T, log_T,
                               sqrt_Is, inverse_T_sqrt_Is, log_T_sqrt_Is, 
                               Is, inverse_T_Is, log_T_Is, inverse_T_sqrt_Is³,
                               inverse_T_Is², log_S, 
                               pressure_correction)
    end
end

@inline (c::KS)(T::FT, S, Is = c.ionic_strength(S); P = nothing) where FT = 
    c.pressure_correction(T, P) *
    exp(c.constant
        + c.inverse_T / T
        + c.log_T * log(T)
        + (c.sqrt_Is + c.inverse_T_sqrt_Is / T + c.log_T_sqrt_Is * log(T)) * √Is
        + (c.Is + c.inverse_T_Is / T + c.log_T_Is * log(T)) * Is
        + c.inverse_T_sqrt_Is³ * Is^convert(FT, 1.5) / T
        + c.inverse_T_Is² * Is^2 / T
        + log(1 + c.log_S * S))

summary(::IO, ::KS) = string("Bisulfate dissociation constant")
show(io::IO, c::KS) = print(io, "Bisulfate dissociation constant\n",
    "    ln(kˢ/k°) = $(c.constant)
                + $(c.inverse_T) / T
                + $(c.log_T) log(T)
                + ($(c.sqrt_Is) + $(c.inverse_T_sqrt_Is) / T + $(c.log_T_sqrt_Is) log(T)) √Is
                + ($(c.Is) + $(c.inverse_T_Is) / T + $(c.log_T_Is) log(T)) Is
                + $(c.inverse_T_sqrt_Is³) √Is³ / T
                + $(c.inverse_T_Is²) Is² / T
                + log(1 + $(c.log_S) S)")                                               

"""
    KF(; ionic_strength = IonicStrength(),
          sulfate_constant = KS(; ionic_strength),
          constant = -9.68,
          inverse_T = 874.0,
          sqrt_S =  0.111,
          log_S = 0.0,
          log_S_KS = 0.0,
          pressure_correction = 
             PressureCorrection(; a₀=-9.78, a₁=-0.0090, a₂=-0.000942, b₀=-0.00391, b₁=0.000054))

Parameterisation for hydrogen fluoride dissociation equilibrium constant.

    HF ⇌ F⁻ + H⁺

    Kᶠ = [H⁺][F⁻]/[HF] 

Default values from Perez and Fraga (1987, Mar. Chem., 21, 161–168).
"""
struct KF{IS, KS, FT, PC}
         ionic_strength :: IS 
       sulfate_constant :: KS 
    
               constant :: FT 
              inverse_T :: FT 
                 sqrt_S :: FT 
                  log_S :: FT 
               log_S_KS :: FT 

    pressure_correction :: PC

    function KF(FT = Float64;
                ionic_strength::IS = IonicStrength{FT}(),
                sulfate_constant::KS = KS(FT; ionic_strength),
                
                constant = -9.68,
                inverse_T = 874.0,
                sqrt_S = 0.111,
                log_S = 0.0,
                log_S_KS = 0.0,

                pressure_correction::PC =
                    PressureCorrection(FT; a₀=-9.78, a₁=-0.0090, a₂=-0.000942, b₀=-0.00391, b₁=0.000054)) where {IS, KS, PC}

        return new{IS, KS, FT, PC}(ionic_strength, sulfate_constant, constant, inverse_T, sqrt_S, log_S, log_S_KS, pressure_correction)
    end

end

@inline (c::KF)(T, S, Is = c.ionic_strength(S), KS = c.sulfate_constant(T, S, Is); P = nothing) =
    c.pressure_correction(T, P) *
    exp(c.constant
        + c.inverse_T / T
        + c.sqrt_S * √S
        + log(1 + c.log_S * S)
        + log(1 + c.log_S_KS * S / KS))

summary(::IO, ::KF) = string("Hydrogen fluoride dissociation constant")
show(io::IO, c::KF) = print(io, "Hydrogen fluoride dissociation constant\n",
    "    ln(kᶠ/k°) = $(c.constant)
                + $(c.inverse_T) / T
                + $(c.sqrt_S) √S
                + log(1 + $(c.log_S) S)
                + log(1 + $(c.log_S_KS) S / Kˢ)")   

"""
    KP(constant,
       inverse_T,
       log_T,
       sqrt_S,
       inverse_T_sqrt_S,
       S,
       inverse_T_S,
       pressure_correction)

Generic equilibrium constant parameterisation of the form used by 
Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677) for phosphoric 
acid dissociation.
"""
struct KP{FT, PC}
               constant :: FT
              inverse_T :: FT
                  log_T :: FT
                 sqrt_S :: FT
       inverse_T_sqrt_S :: FT
                      S :: FT
            inverse_T_S :: FT

    pressure_correction :: PC
end

@inline (c::KP)(T, S; P = nothing) = 
    c.pressure_correction(T, P) *
    exp(c.constant
        + c.inverse_T / T
        + c.log_T * log(T)
        + (c.sqrt_S + c.inverse_T_sqrt_S / T) * √S
        + (c.S + c.inverse_T_S / T) * S)

summary(::IO, ::KP) = string("Phosphate dissociation constant")
show(io::IO, c::KP) = print(io, "Phosphate dissociation constant\n",
    "    ln(kᵖⁿ/k°) = $(c.constant)
                + $(c.inverse_T) / T
                + $(c.log_T) log(T)
                + ($(c.sqrt_S) + $(c.inverse_T_sqrt_S) / T) √S
                + ($(c.S) + $(c.inverse_T_S) / T) S")   

"""
    KP1(; constant = 115.525,
          inverse_T = -4576.752,
          log_T = - 18.453,
          sqrt_S = 0.69171,
          inverse_T_sqrt_S = -106.736,
          S = -0.01844,
          inverse_T_S = -0.65643,
          pressure_correction = 
             PressureCorrection(; a₀=-14.51, a₁=0.1211, a₂=-0.000321, b₀=-0.00267, b₁=0.0000427))

Instance of `KP` returning the first phosphocic acid equilibrium constant.

    H₃PO₄ ⇌ H₂PO₄⁻ + H⁺

    Kᵖ¹ = [H⁺][H₂PO₄]/[H₃PO₄] 

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
KP1(FT = Float64; 
    constant = 115.525,
    inverse_T = -4576.752,
    log_T = - 18.453,
    sqrt_S = 0.69171,
    inverse_T_sqrt_S = -106.736,
    S = -0.01844,
    inverse_T_S = -0.65643,
    pressure_correction::PC =
        PressureCorrection(FT; a₀=-14.51, a₁=0.1211, a₂=-0.000321, b₀=-0.00267, b₁=0.0000427)) where PC = 
    KP{FT, PC}(constant,
               inverse_T, 
               log_T, 
               sqrt_S, 
               inverse_T_sqrt_S, 
               S, 
               inverse_T_S, 
               pressure_correction) 

"""
    KP2(; constant = 172.0883,
          inverse_T = -8814.715,
          log_T = -27.927,
          sqrt_S = 1.3566,
          inverse_T_sqrt_S = -160.340,
          S = -0.05778,
          inverse_T_S = 0.37335,
          pressure_correction = 
            PressureCorrection(; a₀=-23.12, a₁=0.1758, a₂=-0.002647, b₀=-0.00515, b₁=0.00009))

Instance of `KP` returning the second phosphocic acid equilibrium constant.

    H₂PO₄⁻ ⇌ HPO₄²⁻ + H⁺

    Kᵖ² = [H⁺][HPO₄²⁻]/[H₂PO₄⁻] 

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
KP2(FT = Float64; 
    constant = 172.0883,
    inverse_T = -8814.715,
    log_T = -27.927,
    sqrt_S = 1.3566,
    inverse_T_sqrt_S = -160.340,
    S = -0.05778,
    inverse_T_S = 0.37335,
    pressure_correction::PC = 
        PressureCorrection(FT; a₀=-23.12, a₁=0.1758, a₂=-0.002647, b₀=-0.00515, b₁=0.00009)) where PC = 
    KP{FT, PC}(constant,
               inverse_T,
               log_T,
               sqrt_S,
               inverse_T_sqrt_S,
               S,
               inverse_T_S,
               pressure_correction)

"""
    KP3(; constant = -18.141,
          inverse_T = -3070.75,
          log_T = 0.0,
          sqrt_S = 2.81197,
          inverse_T_sqrt_S = 17.27039,
          S = -0.09984,
          inverse_T_S = -44.99486,
          pressure_correction = 
            PressureCorrection(; a₀=-26.57, a₁=0.2020, a₂=-0.0030420, b₀=-0.00408, b₁=0.0000714))

Instance of `KP` returning the third phosphocic acid equilibrium constant.

    HPO₄²⁻ ⇌ PO₄ + H⁺

    Kᵖ³ = [H⁺][PO₄³⁻]/[HPO₄⁻] 

Default values from Millero (1995, Geochim. Cosmochim. Acta, 59, 661–677).
"""
KP3(FT = Float64; 
    constant = - 18.141,
    inverse_T = -3070.75,
    log_T = 0.0,
    sqrt_S = 2.81197,
    inverse_T_sqrt_S = 17.27039,
    S = -0.09984,
    inverse_T_S = -44.99486,
    pressure_correction::PC =
        PressureCorrection(FT; a₀=-26.57, a₁=0.2020, a₂=-0.0030420, b₀=-0.00408, b₁=0.0000714)) where PC = 
    KP{FT, PC}(constant,
               inverse_T,
               log_T,
               sqrt_S,
               inverse_T_sqrt_S,
               S,
               inverse_T_S,
               pressure_correction)

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
struct KSi{IS, FT}
     ionic_strength :: IS 
   
           constant :: FT
          inverse_T :: FT 
              log_T :: FT 
            sqrt_Is :: FT
  inverse_T_sqrt_Is :: FT
                 Is :: FT 
       inverse_T_Is :: FT 
                Is² :: FT 
      inverse_T_Is² :: FT
              log_S :: FT

    function KSi(FT = Float64;
                 ionic_strength::IS = IonicStrength{FT}(),
                 constant =  117.385,
                 inverse_T = -8904.2,
                 log_T = -19.334,
                 sqrt_Is = 3.5913,
                 inverse_T_sqrt_Is = -458.79,
                 Is = -1.5998,
                 inverse_T_Is = 188.74,
                 Is² = 0.07871,
                 inverse_T_Is² = -12.1652,
                 log_S = -0.001005) where IS

        return new{IS, FT}(ionic_strength, constant, inverse_T, log_T, sqrt_Is, inverse_T_sqrt_Is, Is,
                           inverse_T_Is, Is², inverse_T_Is², log_S)
    end
end

@inline (c::KSi)(T, S, Is = c.ionic_strength(S); P = nothing) = 
    exp(c.constant
        + c.inverse_T / T
        + c.log_T * log(T)
        + (c.sqrt_Is + c.inverse_T_sqrt_Is / T) * √Is
        + (c.Is + c.inverse_T_Is / T) * Is
        + (c.Is² + c.inverse_T_Is² / T) * Is^2
        + log(1 + c.log_S * S))

summary(::IO, ::KSi) = string("Silicic acid constant")
show(io::IO, c::KSi) = print(io, "Silicic acid constant\n",
    "    ln(kˢⁱ/k°) = $(c.constant)
                + $(c.inverse_T) / T
                + $(c.log_T) log(T)
                + ($(c.sqrt_Is) + $(c.inverse_T_sqrt_Is) / T) √Is
                + ($(c.Is) + $(c.inverse_T_Is) / T) Is
                + ($(c.Is²) + $(c.inverse_T_Is²) / T) Is²
                + log(1 + $(c.log_S) S)")   

"""
    KSP(therm_constant,
        therm_T,
        therm_inverse_T,
        therm_log_T,
        sea_sqrt_S,
        sea_T_sqrt_S,
        sea_inverse_T_sqrt_S,
        sea_S,
        sea_S_sqrt_S³,
        pressure_correction)

Generic CaCO₃ solubility parameterisation of the form given by
Form from Millero, F. J. (2007, Chemical Reviews, 107(2), 308–341).
"""
struct KSP{FT, PC}
         therm_constant :: FT
                therm_T :: FT
        therm_inverse_T :: FT
            therm_log_T :: FT
             sea_sqrt_S :: FT
           sea_T_sqrt_S :: FT
   sea_inverse_T_sqrt_S :: FT
                  sea_S :: FT
          sea_S_sqrt_S³ :: FT

    pressure_correction :: PC
end

@inline function (c::KSP)(T::FT, S; P = nothing) where FT

    pressure_correction = c.pressure_correction(T, P)

    # The Mucci (1983) / Millero (2007) formula gives log₁₀(KSP), so the
    # final step is 10^x, not exp(x).
    log10KSP = c.therm_constant + c.therm_T * T + c.therm_inverse_T / T + c.therm_log_T * log10(T)
    log10KSP += ((c.sea_sqrt_S + c.sea_T_sqrt_S * T + c.sea_inverse_T_sqrt_S / T) * √S
                  + c.sea_S * S
                  + c.sea_S_sqrt_S³ * S^convert(FT, 1.5))

    return pressure_correction * convert(FT, 10)^log10KSP
end

summary(::IO, ::KSP) = string("Calcite solubility")
show(io::IO, c::KSP) = print(io, "Calcite solubility\n",
    "    log₁₀(kₛₚ) = $(c.therm_constant) + $(c.therm_T) T + $(c.therm_inverse_T) / T + $(c.therm_log_T) log₁₀(T)\n",
    "    log₁₀(kₛₚˢ) = log₁₀(kₛₚ) + ($(c.sea_sqrt_S) + $(c.sea_T_sqrt_S) T + $(c.sea_inverse_T_sqrt_S) / T) √S\n",
    "                + $(c.sea_S) S + $(c.sea_S_sqrt_S³) √S³")   

"""
    KSP_calcite(; therm_constant = -171.9065,
                  therm_T = -0.077993,
                  therm_inverse_T = 2839.319,
                  therm_log_T = 71.595,
                  sea_sqrt_S = -0.77712,
                  sea_T_sqrt_S = 0.0028426,
                  sea_inverse_T_sqrt_S = 178.34,
                  sea_S = -0.07711,
                  sea_S_sqrt_S³ = 0.0041249,
                  pressure_correction =
                      PressureCorrection(; a₀=-48.76, a₁=0.5304, a₂=-0.0, b₀=-0.01176, b₁=0.0003692))

Instance of `KSP` returning calcite solubility.

Default values from Millero, F. J. (2007, Chemical Reviews, 107(2), 308–341).
"""
KSP_calcite(FT = Float64; 
            therm_constant = -171.9065,
            therm_T = -0.077993,
            therm_inverse_T = 2839.319,
            therm_log_T = 71.595,
            sea_sqrt_S = -0.77712,
            sea_T_sqrt_S = 0.0028426,
            sea_inverse_T_sqrt_S = 178.34,
            sea_S = -0.07711,
            sea_S_sqrt_S³ = 0.0041249,
            pressure_correction::PC =
                  PressureCorrection(FT; a₀=-48.76, a₁=0.5304, a₂=-0.0, b₀=-0.01176, b₁=0.0003692)) where PC =
    KSP{FT, PC}(therm_constant,
                therm_T,
                therm_inverse_T,
                therm_log_T,
                sea_sqrt_S,
                sea_T_sqrt_S,
                sea_inverse_T_sqrt_S,
                sea_S,
                sea_S_sqrt_S³,
                pressure_correction)

"""
KSP_aragonite(; therm_constant = -171.945,
                therm_T = -0.077993,
                therm_inverse_T = 2903.293,
                therm_log_T = 71.595,
                sea_sqrt_S = -0.068393,
                sea_T_sqrt_S = 0.0017276,
                sea_inverse_T_sqrt_S = 88.135,
                sea_S = -0.10018,
                sea_S_sqrt_S³ = 0.0059415,
                pressure_correction =
                    PressureCorrection(; a₀=-45.96, a₁=0.5304, a₂=-0.0, b₀=-0.01176, b₁=0.0003692))

Instance of `KSP` returning calcite solubility.

Default values from Millero, F. J. (2007, Chemical Reviews, 107(2), 308–341).
"""
KSP_aragonite(FT = Float64; 
              therm_constant = -171.945,
              therm_T = -0.077993,
              therm_inverse_T = 2903.293,
              therm_log_T = 71.595,
              sea_sqrt_S = -0.068393,
              sea_T_sqrt_S = 0.0017276,
              sea_inverse_T_sqrt_S = 88.135,
              sea_S = -0.10018,
              sea_S_sqrt_S³ = 0.0059415,
              pressure_correction::PC =
                    PressureCorrection(FT; a₀=-45.96, a₁=0.5304, a₂=-0.0, b₀=-0.01176, b₁=0.0003692)) where PC =
    KSP{FT, PC}(therm_constant,
                therm_T,
                therm_inverse_T,
                therm_log_T,
                sea_sqrt_S,
                sea_T_sqrt_S,
                sea_inverse_T_sqrt_S,
                sea_S,
                sea_S_sqrt_S³,
                pressure_correction)
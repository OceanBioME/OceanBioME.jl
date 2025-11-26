using OceanBioME: DampedNewtonRaphsonSolver
using OceanBioME.Models: teos10_polynomial_approximation

struct CarbonChemistry{P0, PC, PB, PS, PF, PP, PSi, PW, IS, PKS, PRho, FV, CV, SO}
          ionic_strength :: IS
              solubility :: P0
           carbonic_acid :: PC
              boric_acid :: PB
                   water :: PW
                 sulfate :: PS
                fluoride :: PF
         phosphoric_acid :: PP
            silicic_acid :: PSi
      calcite_solubility :: PKS
        density_function :: PRho   
first_virial_coefficient :: FV
cross_virial_coefficient :: CV             
                  solver :: SO     
end

"""
    CarbonChemistry(FT = Float64; 
                    ionic_strength = IonicStrength(),
                    solubility = K0(),
                    carbonic_acid = (K1 = K1(), K2 = K2()),
                    boric_acid = KB(),
                    water = KW(),
                    sulfate = KS(; ionic_strength),
                    fluoride = KF(; ionic_strength),
                    phosphoric_acid = (KP1 = KP1(), KP2 = KP2(), KP3 = KP3()),
                    silicic_acid = KSi(; ionic_strength),
                    calcite_solubility = KSP_calcite(),
                    density_function = teos10_polynomial_approximation,
                    first_virial_coefficient = PolynomialVirialCoefficientForCarbonDioxide(),
                    cross_viral_coefficient = CrossVirialCoefficientForCarbonDioxide(),
                    solver = NewtonRaphsonSolver())

Carbon chemistry model capable of solving for sea water pCO₂ from DIC and 
total alkalinity or DIC and pH. 

Default form from Dickson, A.G., Sabine, C.L. and  Christian, J.R. (2007), 
Guide to Best Practices for Ocean CO 2 Measurements. PICES Special Publication 3, 191 pp.

See each parameters documentation for origional sources.

Example
=======

```jldoctest
julia> using OceanBioME

julia> carbon_chemistry = CarbonChemistry()
`CarbonChemistry` model which solves for pCO₂ and pH

julia> pCO₂ = carbon_chemistry(; DIC = 2000.0, Alk = 2000.0, T = 10.0, S = 35.0)
1308.1474527899106

julia> pH = carbon_chemistry(; DIC = 2000.0, Alk = 2000.0, T = 10.0, S = 35.0, output = Val(:pHᶠ))
7.502532746463654

julia> pCO₂_higher_pH = carbon_chemistry(; DIC = 2000.0, T = 10.0, S = 35.0, pH = 7.5)
1315.7136384737507

```
"""
function CarbonChemistry(FT = Float64; 
                         ionic_strength = IonicStrength{FT}(),
                         solubility = K0{FT}(),
                         carbonic_acid = (K1 = K1(FT), K2 = K2(FT)),
                         boric_acid = KB(FT),
                         water = KW(FT),
                         sulfate = KS(FT; ionic_strength),
                         fluoride = KF(FT; ionic_strength),
                         phosphoric_acid = (KP1 = KP1(FT), KP2 = KP2(FT), KP3 = KP3(FT)),
                         silicic_acid = KSi(FT; ionic_strength),
                         calcite_solubility = KSP_calcite(FT),
                         density_function = teos10_polynomial_approximation, # the denisity function *is* going to cause type instability but I can't see a way to fix it
                         first_virial_coefficient = PolynomialVirialCoefficientForCarbonDioxide{FT}(),
                         cross_viral_coefficient = CrossVirialCoefficientForCarbonDioxide{FT}(),

                         solver = DampedNewtonRaphsonSolver{FT, Int, @NamedTuple{lower::FT, upper::FT}}(bounds = (lower = 0, upper = nothing)))

    return CarbonChemistry(ionic_strength, solubility, carbonic_acid, boric_acid, water,
                           sulfate, fluoride, phosphoric_acid, silicic_acid, calcite_solubility, density_function,
                           first_virial_coefficient, cross_viral_coefficient,
                            solver)
end

"""
    (p::CarbonChemistry)(; DIC, T, S, Alk = 0, pH = nothing,
                           output = Val(:fCO₂),
                           boron = 0.000232 / 10.811 * S / 1.80655,
                           sulfate = 0.14 / 96.06 * S / 1.80655,
                           fluoride = 0.000067 / 18.9984 * S / 1.80655,
                           silicate = 0,
                           phosphate = 0,
                           initial_pH_guess = 8)

Calculates `fCO₂` in sea water with `DIC`, `Alk`alinity, `T`emperature, and `S`alinity
unless `pH` is specified, in which case intermediate computation of `pH` is skipped and
`pCO₂` is calculated from the `DIC`, `T`, `S` and `pH`.

`DIC` is expected in mmol C/m³, `Alk` meq/m³, `T` in °C, and S in PSU.

When pH is specified the free pH (i.e. -log[H⁺]) is expected.

Alternativly pCO₂, and free, total, or sea water pH may be returned by setting `output`
to Val(:pCO₂), Val(:pHᶠ), Val(:pHᵗ), or Val(:pHˢ), which will return `X` in `Val(:X)`
instead of fCO₂.
"""
@inline function (p::CarbonChemistry)(; DIC::FT, T, S, Alk = zero(DIC), pH = nothing,
                                        P = nothing, # bars (???)
                                        lon = zero(DIC),
                                        lat = zero(DIC),
                                        output = Val(:fCO₂),
                                        boron = convert(typeof(DIC), 0.000232 / 10.811 * S / 1.80655),
                                        sulfate = convert(typeof(DIC), 0.14 / 96.06 * S / 1.80655),
                                        fluoride = convert(typeof(DIC), 0.000067 / 18.9984 * S / 1.80655),
                                        silicate = zero(DIC),
                                        phosphate = zero(DIC),
                                        initial_pH_guess = convert(typeof(DIC), 8)) where FT

    ρₒ = p.density_function(T, S, ifelse(isnothing(P), one(DIC), P), lon, lat)

    # Centigrade to kelvin
    T += convert(FT, 273.15)

    # mili-equivalents / m³ to equivalents / kg
    Alk *= convert(FT, 1e-3) / ρₒ

    # mmol / m³ to mol / kg
    DIC       *= convert(FT, 1e-3) / ρₒ
    phosphate *= convert(FT, 1e-3) / ρₒ
    silicate  *= convert(FT, 1e-3) / ρₒ
    
    # ionic strength
    Is = p.ionic_strength(S)

    # compute equilibrium constants
    K1  = p.carbonic_acid.K1(T, S; P)
    K2  = p.carbonic_acid.K2(T, S; P)
    KB  = p.boric_acid(T, S; P)
    KW  = p.water(T, S; P)
    KS  = p.sulfate(T, S, Is; P)
    KF  = p.fluoride(T, S, Is, KS; P)
    KP1 = p.phosphoric_acid.KP1(T, S; P)
    KP2 = p.phosphoric_acid.KP2(T, S; P)
    KP3 = p.phosphoric_acid.KP3(T, S; P)
    KSi = p.silicic_acid(T, S, Is)

    params = (; DIC, Alk, boron, sulfate, fluoride, silicate, phosphate,
                K1, K2, KB, KW, KS, KF, KP1, KP2, KP3, KSi)

    # solve equilibrium for hydrogen ion concentration
    H = solve_for_H(pH, params, initial_pH_guess, p.solver)

    # compute solubility equilibrium constant
    K0 = p.solubility(T, S)

    # compute pCO₂
    CO₂  = DIC * H ^ 2 / (H ^ 2 + K1 * H + K1 * K2) 
    fCO₂ = (CO₂ / K0) * convert(FT, 10 ^ 6) # μatm

    return selected_output(output, fCO₂, H, P, T, S, Is, sulfate, fluoride, p)
end

@inline selected_output(::Val{:fCO₂}, fCO₂, H, P, Tk, S, Is, sulfate, fluoride, p) = fCO₂ # ppm
@inline selected_output(::Val{:pHᶠ}, fCO₂, H, P, Tk, S, Is, sulfate, fluoride, p) = -log10(H) # 

@inline function selected_output(::Val{:pCO₂}, fCO₂::FT, H, P, Tk, S, Is, sulfate, fluoride, p) where FT
    P = ifelse(isnothing(P), one(fCO₂), P)
    P *= convert(FT, 101325) # pascals

    B = p.first_virial_coefficient(Tk)
    δ = p.cross_virial_coefficient(Tk)

    fCO₂ *= convert(FT, 0.09807) # μatm -> Pa

    φ = one(FT)

    xCO₂ = fCO₂ / (φ * P) # Pa to mol/mol

    # Experimentally this converged xCO₂ to machine precision
    for n = 1:3
        φ = exp((B + convert(FT, 2) * (one(FT) - xCO₂)^convert(FT, 2) * δ) * P / (convert(FT, GAS_CONSTANT) * Tk))
        xCO₂ = fCO₂ / (φ * P)
    end

    pCO₂ = fCO₂ / φ # Pa
    pCO₂ /= convert(FT, 0.09807) # μatm or ppmv

    return pCO₂ # ppmv
end

@inline function selected_output(::Val{:pHᵗ}, fCO₂::FT, H, P, Tk, S, Is, sulfate, fluoride, p) where FT
    KS  = p.sulfate(Tk, S, Is; P)
    HSO₄⁻ = sulfate / (1 + KS / H)

    return -log10(H + HSO₄⁻)
end

@inline function selected_output(::Val{:pHˢ}, fCO₂::FT, H, P, Tk, S, Is, sulfate, fluoride, p) where FT
    KS  = p.sulfate(Tk, S, Is; P)
    HSO₄⁻ = sulfate / (1 + KS / H)

    KF  = p.fluoride(Tk, S, Is, KS; P)
    HF = fluoride / (1 + KF / H)

    return -log10(H + HSO₄⁻ + HF)
end

# solves `alkalinity_residual` for pH
solve_for_H(pH::FT, args...) where FT = convert(FT, 10.0) ^ - pH

include("alkalinity_residual.jl")

solve_for_H(::Nothing, params, initial_pH_guess::FT, solver) where FT =
    solver(alkalinity_residual, ∂ₕ_alkalinity_residual, 10^(-initial_pH_guess), params)#

# display
summary(::IO, ::CarbonChemistry) = string("`CarbonChemistry` model")
show(io::IO, ::CarbonChemistry) = print(io, "`CarbonChemistry` model which solves for pCO₂ and pH")

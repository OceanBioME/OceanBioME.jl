"""
    PartiallySolubleGas(; air_concentration, solubility)

Parameterises the available concentration of a gas dissolving in water in the form ``\\alpha C_a``
where ``\alpha`` is the Ostwald solubility coeffieient and ``C_a`` is the concentration in the air.
"""
struct PartiallySolubleGas{AC, S}
    air_concentration :: AC
           solubility :: S

    function PartiallySolubleGas(FT = Float64;
                                 air_concentration,
                                 solubility :: S) where S

        air_concentration = normalise_surface_function(air_concentration; FT)

        AC = typeof(air_concentration)

        return new{AC, S}(air_concentration, solubility)
    end

    # bypasses `normalise_surface_function`, for `Adapt`
    PartiallySolubleGas{AC, S}(air_concentration::AC, solubility::S) where {AC, S} =
        new{AC, S}(air_concentration, solubility)
end

Adapt.adapt_structure(to, gs::PartiallySolubleGas) =
    PartiallySolubleGas{typeof(adapt(to, gs.air_concentration)),
                        typeof(adapt(to, gs.solubility))}(adapt(to, gs.air_concentration),
                                                          adapt(to, gs.solubility))

@inline surface_value(gs::PartiallySolubleGas, i, j, grid, clock, model_fields) = 
    surface_value(gs.air_concentration, i, j, grid, clock) * surface_value(gs.solubility, i, j, grid, clock, model_fields)

"""
    Wanninkhof92Solubility

Parameterises the Ostwald solubility coefficient as given in Wanninkhof, 1992.
"""
struct Wanninkhof92Solubility{FT}
    A1 :: FT
    A2 :: FT
    A3 :: FT
    B1 :: FT
    B2 :: FT
    B3 :: FT
end

function surface_value(sol::Wanninkhof92Solubility, i, j, grid, clock, model_fields)
    FT = eltype(grid)

    Tk = @inbounds model_fields.T[i, j, grid.Nz] + convert(FT, 273.15)
    S = @inbounds model_fields.S[i, j, grid.Nz]

    Tk_100 = Tk / convert(FT, 100)

    β = exp(sol.A1 + sol.A2 / Tk_100 + sol.A3 * log(Tk_100) + S * (sol.B1 + sol.B2 * Tk_100 + sol.B2 * Tk_100^convert(FT, 2)))

    return β / Tk
end

OxygenSolubility(FT = Float64; A1 = -58.3877, A2 = 85.8079, A3 = 23.8439, B1 = -0.034892, B2 = 0.015578, B3 = -0.0019387) =
    Wanninkhof92Solubility{FT}(A1, A2, A3, B1, B2, B3)

struct MolPerKgPerAtmToMMolPerCubicMPerMicroAtm{SO, DE}
    solubility :: SO
       density :: DE
end

Adapt.adapt_structure(to, k::MolPerKgPerAtmToMMolPerCubicMPerMicroAtm) = 
    MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(adapt(to, k.solubility),
                                             adapt(to, k.density))

@inline (dss::MolPerKgPerAtmToMMolPerCubicMPerMicroAtm)(T::FT, S) where FT = dss.solubility(T + FT(273.15), S) * dss.density(T, S) / FT(10^3)


"""
    GarciaGordonOxygenSaturation

The saturation concentration of oxygen in sea water (mmol O₂ / m³) after Garcia and Gordon
(1992), Limnology and Oceanography, page 1310, equation (8), scaled by the atmospheric
pressure. Note that the `A₃Tₛ²` term printed in the paper is an error and is not included
(as in the reference implementation).

The fit is written in terms of the scaled temperature
``T_s = \\ln\\left[(T_0 + T_{ref} - T)/(T_0 + T)\\right]`` (with ``T_0`` the freezing point
of fresh water in Kelvin and ``T_{ref}`` a reference temperature in °C), as

```math
[O_2]_{sat} = \\frac{p_{atm}}{V_m}\\exp\\left[A(T_s) + S\\left(B(T_s) + Sc_0\\right)\\right],
```

where ``A`` and ``B`` are polynomials in ``T_s``, ``V_m`` is the molar volume of oxygen
(which converts the ml/l of the original fit to mmol/m³), and ``p_{atm}`` is the
atmospheric pressure in atmospheres.

The fit is valid for ``T_{freezing} \\leq T \\leq 40°C`` and ``0 \\leq S \\leq 42``; it is
not clamped to that range. Its check value is 282.015 mmol/m³ at ``T = 10°C``, ``S = 35``,
and ``p_{atm} = 1``.

This is an alternative to a [`PartiallySolubleGas`](@ref) for the `air_concentration` of an
[`OxygenGasExchangeBoundaryCondition`](@ref).
"""
struct GarciaGordonOxygenSaturation{FT, TP, SP, AP}
          temperature_polynomial :: TP  # MARBL a_0 … a_5
             salinity_polynomial :: SP  # MARBL b_0 … b_3
    salinity_squared_coefficient :: FT  # MARBL c_0
           temperature_reference :: FT  # °C, part of the definition of the fit
             oxygen_molar_volume :: FT  # m³/mol, converts the fit's ml/l to mmol/m³
            atmospheric_pressure :: AP  # atm; a scalar, function, or `Field`
end

"""
    GarciaGordonOxygenSaturation(FT = Float64;
                                 temperature_coefficients = (2.00907, 3.22014, 4.05010, 4.94457, -2.56847e-1, 3.88767),
                                 salinity_coefficients = (-6.24523e-3, -7.37614e-3, -1.03410e-2, -8.17083e-3),
                                 salinity_squared_coefficient = -4.88682e-7,
                                 temperature_reference = 25.0,       # °C
                                 oxygen_molar_volume = 0.0223916,    # m³/mol
                                 atmospheric_pressure = 1)           # atm

Returns the Garcia and Gordon (1992) oxygen saturation concentration in mmol O₂ / m³.

Keyword Arguments
=================

- `temperature_coefficients`: the six coefficients of the polynomial in the scaled
  temperature ``T_s``
- `salinity_coefficients`: the four coefficients of the polynomial in ``T_s`` multiplying
  the salinity
- `salinity_squared_coefficient`: the coefficient of the salinity squared term
- `temperature_reference`: the reference temperature of the scaled temperature (°C)
- `oxygen_molar_volume`: the molar volume of oxygen (m³/mol), which converts the ml/l of
  the original fit to mmol/m³
- `atmospheric_pressure`: the atmospheric pressure (atm), which may be a number, a function
  of the form `(x, y, t)`, or a `Field`; the default of 1 leaves the saturation at the
  1 atm value of the fit

See also [`PartiallySolubleGas`](@ref) and [`Wanninkhof92Solubility`](@ref).
"""
function GarciaGordonOxygenSaturation(FT = Float64;
                                      temperature_coefficients = (2.00907, 3.22014, 4.05010, 4.94457, -2.56847e-1, 3.88767),
                                      salinity_coefficients = (-6.24523e-3, -7.37614e-3, -1.03410e-2, -8.17083e-3),
                                      salinity_squared_coefficient = -4.88682e-7,
                                      temperature_reference = 25.0,    # °C
                                      oxygen_molar_volume = 0.0223916, # m³/mol
                                      atmospheric_pressure = 1)        # atm

    temperature_polynomial = PolynomialParameterisation{5}(FT; coefficients = temperature_coefficients)
    salinity_polynomial    = PolynomialParameterisation{3}(FT; coefficients = salinity_coefficients)

    atmospheric_pressure = normalise_surface_function(atmospheric_pressure; FT)

    TP = typeof(temperature_polynomial)
    SP = typeof(salinity_polynomial)
    AP = typeof(atmospheric_pressure)

    return GarciaGordonOxygenSaturation{FT, TP, SP, AP}(temperature_polynomial,
                                                        salinity_polynomial,
                                                        convert(FT, salinity_squared_coefficient),
                                                        convert(FT, temperature_reference),
                                                        convert(FT, oxygen_molar_volume),
                                                        atmospheric_pressure)
end

@inline function surface_value(sat::GarciaGordonOxygenSaturation, i, j, grid, clock, model_fields)
    FT = eltype(grid)

    T = @inbounds model_fields.T[i, j, grid.Nz]
    S = @inbounds model_fields.S[i, j, grid.Nz]

    T₀ = convert(FT, 273.15)

    Tₛ = log((T₀ + sat.temperature_reference - T) / (T₀ + T))

    A = sat.temperature_polynomial(Tₛ)
    B = sat.salinity_polynomial(Tₛ)

    O₂sat = exp(A + S * (B + S * sat.salinity_squared_coefficient)) / sat.oxygen_molar_volume

    return O₂sat * surface_value(sat.atmospheric_pressure, i, j, grid, clock)
end

Adapt.adapt_structure(to, sat::GarciaGordonOxygenSaturation{FT}) where FT =
    GarciaGordonOxygenSaturation{FT,
                                 typeof(adapt(to, sat.temperature_polynomial)),
                                 typeof(adapt(to, sat.salinity_polynomial)),
                                 typeof(adapt(to, sat.atmospheric_pressure))}(adapt(to, sat.temperature_polynomial),
                                                                              adapt(to, sat.salinity_polynomial),
                                                                              sat.salinity_squared_coefficient,
                                                                              sat.temperature_reference,
                                                                              sat.oxygen_molar_volume,
                                                                              adapt(to, sat.atmospheric_pressure))

summary(::GarciaGordonOxygenSaturation) = "GarciaGordonOxygenSaturation"

show(io::IO, sat::GarciaGordonOxygenSaturation) =
    println(io, summary(sat), "\n",
                "    [O₂]sat(T, S) = pₐₜₘ exp(A(Tₛ) + S(B(Tₛ) + Sc₀)) / Vₘ, Tₛ = ln((T₀ + $(sat.temperature_reference) - T)/(T₀ + T)),\n",
                "    A: $(sat.temperature_polynomial.coefficients),\n",
                "    B: $(sat.salinity_polynomial.coefficients),\n",
                "    c₀ = $(sat.salinity_squared_coefficient), Vₘ = $(sat.oxygen_molar_volume) m³/mol, pₐₜₘ = $(sat.atmospheric_pressure) atm")

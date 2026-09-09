"""
    CarbonDioxideConcentration(FT = Float64;
                               carbon_chemistry::CC,
                               DIC = :DIC,
                               Alk = :Alk)

The water-side carbon dioxide concentration seen by an air-sea gas exchange: the aqueous
carbon dioxide concentration, `[CO₂(aq)]` in mmol / m³, computed by the `carbon_chemistry`
model from the model's dissolved inorganic carbon and alkalinity.

`DIC` and `Alk` specify the tracer names of the DIC and alkalinity in the model.

The exchange is a difference of concentrations: the solubility that turns the air-side mole
fraction into a concentration lives on the `air_concentration`
(see [`CarbonDioxideAirConcentration`](@ref)), and the transfer velocity is a bare piston
velocity.

Follows Dickson, A.G., Sabine, C.L. and Christian, J.R. (2007), Guide to Best Practices for
Ocean CO 2 Measurements. PICES Special Publication 3, 191 pp.

No atmospheric pressure enters the water-side concentration (see
[`CarbonDioxideAirConcentration`](@ref), which carries the only atmospheric pressure on this
path).
"""
struct CarbonDioxideConcentration{DIC, Alk, CC<:CarbonChemistry}
    carbon_chemistry :: CC
end

CarbonDioxideConcentration(FT = Float64;
                           carbon_chemistry::CC = CarbonChemistry(FT),
                           DIC = :DIC,
                           Alk = :Alk) where CC =
    CarbonDioxideConcentration{DIC, Alk, CC}(carbon_chemistry)

summary(::CarbonDioxideConcentration{DIC, Alk, CC}) where {DIC, Alk, CC} =
    "`CarbonChemistry` derived aqueous carbon dioxide concentration ([CO₂(aq)], mmol/m³) {$DIC, $Alk, $(nameof(CC))}"

show(io::IO, ccc::CarbonDioxideConcentration{DIC, Alk}) where {DIC, Alk} =
    println(io, summary(ccc), "\n",
            "    Solves the $(nameof(typeof(ccc.carbon_chemistry))) based on $DIC and $Alk")

@inline function surface_value(cc::CarbonDioxideConcentration{DIC_name, Alk_name}, i, j, grid, clock, model_fields) where {DIC_name, Alk_name}
    DIC = @inbounds model_fields[DIC_name][i, j, grid.Nz] # this is a compile time inference so is fine on GPU
    Alk = @inbounds model_fields[Alk_name][i, j, grid.Nz]

    T = @inbounds model_fields.T[i, j, grid.Nz]
    S = @inbounds model_fields.S[i, j, grid.Nz]

    silicate  = silicate_concentration(grid, i, j, grid.Nz, model_fields)
    phosphate = phosphate_concentration(grid, i, j, grid.Nz, model_fields)

    return cc.carbon_chemistry(; DIC, Alk, T, S, silicate, phosphate, output = Val(:CO₂))
end

"""
    CarbonDioxideAirConcentration

The air-side carbon dioxide concentration seen by an air-sea gas exchange, given by Dalton's
law as the product of the dry-air mole fraction and the total atmospheric pressure, and
converted into the units of the water-side concentration by a `solubility`,

```math
x_{CO_2} p_{atm} f_f(T, S) \\rho / 10^3,
```

The default `solubility` is the Weiss and Price (1980) [`FF`](@ref) fit converted to
mmol / m³ per μatm, which pairs with a [`CarbonDioxideConcentration`](@ref) water side. Note
that this is the solubility of a *dry air mole fraction*
(``f_f = K_0 (1 - p_{H_2O}) \\gamma``), and so already carries the water vapour and
non-ideality corrections; it is not [`K0`](@ref).

A `solubility` of `nothing` leaves the air concentration as a mole fraction in ppmv.

This drops into the `air_concentration` of a
[`CarbonDioxideGasExchangeBoundaryCondition`](@ref) in place of a bare number, and with the
default pressure of 1 atm it leaves the mole fraction unscaled by pressure.

Note that this is the *only* atmospheric pressure on the carbon dioxide exchange path. The
water-side [`CarbonDioxideConcentration`](@ref) carries no atmospheric pressure, and the
hydrostatic pressure of the `CarbonChemistry` fugacity coefficient is a distinct quantity
which is not this Dalton factor.
"""
struct CarbonDioxideAirConcentration{MF, AP, SO}
           mole_fraction :: MF # ppmv (≡ μatm at 1 atm)
    atmospheric_pressure :: AP # atm
              solubility :: SO # mmol/m³ per μatm, or `nothing` to stay in ppmv
end

"""
    CarbonDioxideAirConcentration(FT = Float64;
                                  mole_fraction = 413,       # ppmv
                                  atmospheric_pressure = 1,  # atm
                                  solubility = nothing)

Returns the air-side carbon dioxide concentration ``x_{CO_2} p_{atm}`` (in ppmv ≡ μatm), or
``x_{CO_2} p_{atm} f_f \\rho / 10^3`` (in mmol / m³) when a `solubility` is given.

Keyword Arguments
=================

- `mole_fraction`: the dry-air mole fraction of carbon dioxide (ppmv), which may be a
  number, a function of the form `(x, y, t)`, or a `Field`
- `atmospheric_pressure`: the total atmospheric pressure (atm), which may be a number, a
  function of the form `(x, y, t)`, or a `Field`; the default of 1 leaves the mole fraction
  unchanged. Note that the units are atmospheres, not pascals; this is checked (with a
  warning) only when a number is given, since the value of a function or `Field` is not
  known at construction time
- `solubility`: a function of `(T, S)` returning the conversion from a partial pressure in
  μatm to a concentration in mmol / m³, or `nothing` (the default) to leave the air
  concentration as a mole fraction in ppmv. `CarbonDioxideGasExchangeBoundaryCondition`
  supplies `MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(FF{FT}(), density_function)` here

See also [`CarbonDioxideConcentration`](@ref) and
[`CarbonDioxideGasExchangeBoundaryCondition`](@ref).
"""
function CarbonDioxideAirConcentration(FT = Float64;
                                       mole_fraction = 413,      # ppmv
                                       atmospheric_pressure = 1, # atm
                                       solubility = nothing)

    mole_fraction = normalise_surface_function(mole_fraction; FT)
    atmospheric_pressure = normalise_surface_function(atmospheric_pressure; FT)

    MF = typeof(mole_fraction)
    AP = typeof(atmospheric_pressure)
    SO = typeof(solubility)

    return CarbonDioxideAirConcentration{MF, AP, SO}(mole_fraction, atmospheric_pressure, solubility)
end

@inline surface_value(ac::CarbonDioxideAirConcentration, i, j, grid, clock, model_fields) =
    surface_value(ac.mole_fraction, i, j, grid, clock) *
    surface_value(ac.atmospheric_pressure, i, j, grid, clock) *
    air_solubility(ac.solubility, i, j, grid, model_fields)

# no solubility: the air concentration stays a mole fraction in ppmv
@inline air_solubility(::Nothing, i, j, grid, model_fields) = one(eltype(grid))

@inline air_solubility(solubility, i, j, grid, model_fields) =
    solubility(@inbounds(model_fields.T[i, j, grid.Nz]),
               @inbounds(model_fields.S[i, j, grid.Nz]))

Adapt.adapt_structure(to, ac::CarbonDioxideAirConcentration) =
    CarbonDioxideAirConcentration{typeof(adapt(to, ac.mole_fraction)),
                                  typeof(adapt(to, ac.atmospheric_pressure)),
                                  typeof(adapt(to, ac.solubility))}(adapt(to, ac.mole_fraction),
                                                                    adapt(to, ac.atmospheric_pressure),
                                                                    adapt(to, ac.solubility))

summary(::CarbonDioxideAirConcentration{<:Any, <:Any, Nothing}) = "Dalton's law `CarbonDioxideAirConcentration` (xCO₂ pₐₜₘ, ppmv)"
summary(::CarbonDioxideAirConcentration) = "Dalton's law `CarbonDioxideAirConcentration` (xCO₂ pₐₜₘ ff ρ / 10³, mmol/m³)"

show(io::IO, ac::CarbonDioxideAirConcentration) =
    println(io, summary(ac), "\n",
                "    xCO₂ = $(ac.mole_fraction) ppmv,\n",
                "    pₐₜₘ = $(ac.atmospheric_pressure) atm,\n",
                "    solubility = $(ac.solubility)")

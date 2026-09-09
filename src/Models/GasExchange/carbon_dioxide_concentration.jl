"""
    CarbonDioxideConcentration(FT = Float64;
                               carbon_chemistry::CC,
                               DIC = :DIC,
                               Alk = :Alk,
                               output = Val(:CO₂))

The water-side carbon dioxide concentration seen by an air-sea gas exchange, computed by the
`carbon_chemistry` model from the model's dissolved inorganic carbon and alkalinity.

`DIC` and `Alk` specify the tracer names of the DIC and alkalinity in the model.

`output` selects which `carbon_chemistry` quantity is returned, and thereby the *basis* on
which the exchange is computed:

- `Val(:CO₂)` (the default) returns the aqueous carbon dioxide concentration, `[CO₂(aq)]` in
  mmol / m³, so the exchange is a difference of *concentrations*. The solubility which turns
  the air-side mole fraction into a concentration then lives on the `air_concentration`
  (see [`CarbonDioxideAirConcentration`](@ref)), and the transfer velocity is a bare piston
  velocity.
- `Val(:pCO₂)` returns the partial pressure of carbon dioxide in μatm, so the exchange is a
  difference of *partial pressures*, and the solubility instead lives on the
  `transfer_velocity`. This is the legacy basis: it compares a dry air mole fraction against
  a moist air partial pressure, and it passes the water side through the fugacity
  coefficient, which does not cancel. It is retained only to reproduce older results.

[`CarbonDioxideGasExchangeBoundaryCondition`](@ref) chooses matching defaults for the
`transfer_velocity` and `air_concentration` from this choice, so the two bases cannot be
mixed by accident.

Follows Dickson, A.G., Sabine, C.L. and Christian, J.R. (2007), Guide to Best Practices for
Ocean CO 2 Measurements. PICES Special Publication 3, 191 pp.

The fugacity coefficient, and hence the virial coefficients that set it, belong to the
`carbon_chemistry`; and no atmospheric pressure enters the water-side concentration (see
[`CarbonDioxideAirConcentration`](@ref), which carries the only atmospheric pressure on this
path).
"""
struct CarbonDioxideConcentration{DIC, Alk, CC<:CarbonChemistry, Out}
    carbon_chemistry :: CC 
end

CarbonDioxideConcentration(FT = Float64;
                           carbon_chemistry::CC = CarbonChemistry(FT),
                           DIC = :DIC,
                           Alk = :Alk,
                           output::Out = Val(:CO₂)) where {CC, Out} = 
    CarbonDioxideConcentration{DIC, Alk, CC, Out}(carbon_chemistry)

carbon_dioxide_output_name(::Val{:CO₂})  = "aqueous carbon dioxide concentration ([CO₂(aq)], mmol/m³)"
carbon_dioxide_output_name(::Val{:pCO₂}) = "partial pressure of CO₂ (pCO₂, μatm)"
carbon_dioxide_output_name(::Val{name}) where name = string(name)

summary(::CarbonDioxideConcentration{DIC, Alk, CC, Out}) where {DIC, Alk, CC, Out} = 
    "`CarbonChemistry` derived $(carbon_dioxide_output_name(Out())) {$DIC, $Alk, $(nameof(CC))}"

show(io::IO, ccc::CarbonDioxideConcentration{DIC, Alk, CC, Out}) where {DIC, Alk, CC, Out} = 
    println(io, summary(ccc), "\n",
            "    Solves the $(nameof(CC)) based on $DIC and $Alk")

@inline function surface_value(cc::CarbonDioxideConcentration{DIC_name, Alk_name, <:Any, Out}, i, j, grid, clock, model_fields) where {DIC_name, Alk_name, Out}
    DIC = @inbounds model_fields[DIC_name][i, j, grid.Nz] # this is a compile time inference so is fine on GPU
    Alk = @inbounds model_fields[Alk_name][i, j, grid.Nz]

    T = @inbounds model_fields.T[i, j, grid.Nz]
    S = @inbounds model_fields.S[i, j, grid.Nz]

    silicate  = silicate_concentration(grid, i, j, grid.Nz, model_fields)
    phosphate = phosphate_concentration(grid, i, j, grid.Nz, model_fields)

    # mmol/m³ for the default `Val(:CO₂)`, μatm for the legacy `Val(:pCO₂)`
    return cc.carbon_chemistry(; DIC, Alk, T, S, silicate, phosphate, output = Out())
end

"""
    CarbonDioxideAirConcentration

The air-side carbon dioxide concentration seen by an air-sea gas exchange, given by Dalton's
law as the product of the dry-air mole fraction and the total atmospheric pressure, and
converted into the units of the water-side concentration by a `solubility`,

```math
x_{CO_2} p_{atm} f_f(T, S) \\rho / 10^3,
```

matching the way the reference implementation forms its air-side term (the atmospheric
pressure multiplies the air-side concentration, linearly and once, and enters neither the
water-side concentration nor the transfer velocity; and the solubility sits here rather than
on the piston velocity).

The default `solubility` is the Weiss and Price (1980) [`FF`](@ref) fit converted to
mmol / m³ per μatm, which pairs with a water-side
[`CarbonDioxideConcentration`](@ref) returning `Val(:CO₂)`. Note that this is the solubility
of a *dry air mole fraction* (``f_f = K_0 (1 - p_{H_2O}) \\gamma``), and so already carries
the water vapour and non-ideality corrections; it is not [`K0`](@ref).

A `solubility` of `nothing` leaves the air concentration as a mole fraction in ppmv, which
pairs with the legacy `Val(:pCO₂)` water side (there the conversion is carried by the
`transfer_velocity` instead).

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
  supplies `MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(FF{FT}(), density_function)` here when
  the water side is on the concentration basis

See also [`CarbonDioxideConcentration`](@ref) and
[`CarbonDioxideGasExchangeBoundaryCondition`](@ref).
"""
function CarbonDioxideAirConcentration(FT = Float64;
                                       mole_fraction = 413,      # ppmv
                                       atmospheric_pressure = 1, # atm
                                       solubility = nothing)

    if atmospheric_pressure isa Number && atmospheric_pressure > 10
        @warn "The `atmospheric_pressure` $(atmospheric_pressure) is very large, are you sure it is in atmospheres (and not, for example, pascals)?"
    end

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

# no solubility: the air concentration stays a mole fraction in ppmv (the legacy basis)
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

"""
    CarbonDioxideConcentration(FT = Float64;
                               carbon_chemistry::CC,
                               DIC = :DIC,
                               Alk = :Alk)

Converts fCO₂ to partial pressure as per Dickson, A.G., Sabine, C.L. and  Christian, J.R. (2007), 
Guide to Best Practices for Ocean CO 2 Measurements. PICES Special Publication 3, 191 pp.

`DIC` and `Alk` specify the tracer names of the DIC and alkalinity in the model.

The fugacity coefficient, and hence the virial coefficients that set it, belong to the
`carbon_chemistry`; and no atmospheric pressure enters the water-side concentration (see
[`CarbonDioxideAirConcentration`](@ref), which carries the only atmospheric pressure on this
path).
"""
struct CarbonDioxideConcentration{DIC, Alk, CC<:CarbonChemistry}
    carbon_chemistry :: CC 
end

CarbonDioxideConcentration(FT = Float64;
                           carbon_chemistry::CC = CarbonChemistry(),
                           DIC = :DIC,
                           Alk = :Alk) where CC = 
    CarbonDioxideConcentration{DIC, Alk, CC}(carbon_chemistry)

summary(::CarbonDioxideConcentration{DIC, Alk, CC}) where {DIC, Alk, CC} = 
    "`CarbonChemistry` derived partial pressure of CO₂ (pCO₂) {$DIC, $Alk, $(nameof(CC))}"

show(io::IO, ::CarbonDioxideConcentration{DIC, Alk, CC}) where {DIC, Alk, CC} = 
    println(io, "`CarbonChemistry` derived partial pressure of CO₂ (pCO₂) \n",
            "    Solves the $(nameof(CC)) based on $DIC and $Alk")

@inline function surface_value(cc::CarbonDioxideConcentration{DIC_name, Alk_name}, i, j, grid, clock, model_fields) where {DIC_name, Alk_name}
    DIC = @inbounds model_fields[DIC_name][i, j, grid.Nz] # this is a compile time inference so is fine on GPU
    Alk = @inbounds model_fields[Alk_name][i, j, grid.Nz]

    T = @inbounds model_fields.T[i, j, grid.Nz]
    S = @inbounds model_fields.S[i, j, grid.Nz]

    silicate  = silicate_concentration(grid, i, j, grid.Nz, model_fields)
    phosphate = phosphate_concentration(grid, i, j, grid.Nz, model_fields)

    pCO₂ = cc.carbon_chemistry(; DIC, Alk, T, S, silicate, phosphate, output = Val(:pCO₂))

    return pCO₂ # ppmv
end

"""
    CarbonDioxideAirConcentration

The air-side carbon dioxide concentration seen by an air-sea gas exchange, given by
Dalton's law as the product of the dry-air mole fraction and the total atmospheric
pressure,

```math
x_{CO_2} p_{atm},
```

matching the way the reference implementation forms its air-side term (the atmospheric
pressure multiplies the air-side concentration, linearly and once, and enters neither the
water-side concentration nor the transfer velocity).

This drops into the `air_concentration` of a
[`CarbonDioxideGasExchangeBoundaryCondition`](@ref) in place of a bare number, and with the
default pressure of 1 atm it returns the mole fraction unchanged.

Note that this is the *only* atmospheric pressure on the carbon dioxide exchange path. The
water-side [`CarbonDioxideConcentration`](@ref) carries no atmospheric pressure, and the
hydrostatic pressure of the `CarbonChemistry` fugacity coefficient is a distinct quantity
which is not this Dalton factor.
"""
struct CarbonDioxideAirConcentration{MF, AP}
           mole_fraction :: MF # ppmv (≡ μatm at 1 atm)
    atmospheric_pressure :: AP # atm
end

"""
    CarbonDioxideAirConcentration(FT = Float64;
                                  mole_fraction = 413,      # ppmv
                                  atmospheric_pressure = 1) # atm

Returns the air-side carbon dioxide concentration ``x_{CO_2} p_{atm}`` in ppmv (≡ μatm).

Keyword Arguments
=================

- `mole_fraction`: the dry-air mole fraction of carbon dioxide (ppmv), which may be a
  number, a function of the form `(x, y, t)`, or a `Field`
- `atmospheric_pressure`: the total atmospheric pressure (atm), which may be a number, a
  function of the form `(x, y, t)`, or a `Field`; the default of 1 leaves the mole fraction
  unchanged. Note that the units are atmospheres, not pascals; this is checked (with a
  warning) only when a number is given, since the value of a function or `Field` is not
  known at construction time

See also [`CarbonDioxideConcentration`](@ref) and
[`CarbonDioxideGasExchangeBoundaryCondition`](@ref).
"""
function CarbonDioxideAirConcentration(FT = Float64;
                                       mole_fraction = 413,      # ppmv
                                       atmospheric_pressure = 1) # atm

    if atmospheric_pressure isa Number && atmospheric_pressure > 10
        @warn "The `atmospheric_pressure` $(atmospheric_pressure) is very large, are you sure it is in atmospheres (and not, for example, pascals)?"
    end

           mole_fraction = normalise_surface_function(mole_fraction; FT)
    atmospheric_pressure = normalise_surface_function(atmospheric_pressure; FT)

    MF = typeof(mole_fraction)
    AP = typeof(atmospheric_pressure)

    return CarbonDioxideAirConcentration{MF, AP}(mole_fraction, atmospheric_pressure)
end

@inline surface_value(ac::CarbonDioxideAirConcentration, i, j, grid, clock, model_fields) =
    surface_value(ac.mole_fraction, i, j, grid, clock) *
    surface_value(ac.atmospheric_pressure, i, j, grid, clock)

Adapt.adapt_structure(to, ac::CarbonDioxideAirConcentration) =
    CarbonDioxideAirConcentration{typeof(adapt(to, ac.mole_fraction)),
                                  typeof(adapt(to, ac.atmospheric_pressure))}(adapt(to, ac.mole_fraction),
                                                                              adapt(to, ac.atmospheric_pressure))

summary(::CarbonDioxideAirConcentration) = "Dalton's law `CarbonDioxideAirConcentration` (xCO₂ pₐₜₘ)"

show(io::IO, ac::CarbonDioxideAirConcentration) =
    println(io, summary(ac), "\n",
                "    xCO₂ = $(ac.mole_fraction) ppmv,\n",
                "    pₐₜₘ = $(ac.atmospheric_pressure) atm")

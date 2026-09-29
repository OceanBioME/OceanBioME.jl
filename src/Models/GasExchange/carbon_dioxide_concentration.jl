"""
    CarbonDioxideConcentration(FT = Float64;
                               carbon_chemistry::CC,
                               DIC = :DIC,
                               Alk = :Alk,
                               warm_start = false,
                               grid = nothing)

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

With `warm_start = true` the surface free pH from each solve is stored in a
`Field{Center, Center, Nothing}` on `grid` (the `pH` property) and used as the
`initial_pH_guess` of the next, which roughly halves the solver iterations since the surface pH
changes little between time steps. The field starts at zero, and any stored value outside
`0 < pH < 14` (e.g. the first solve, or after a restart) falls back to the default guess of 8.
"""
struct CarbonDioxideConcentration{DIC, Alk, CC<:CarbonChemistry, PH}
    carbon_chemistry :: CC
                  pH :: PH
end

CarbonDioxideConcentration{DIC, Alk}(carbon_chemistry::CC, pH::PH) where {DIC, Alk, CC, PH} =
    CarbonDioxideConcentration{DIC, Alk, CC, PH}(carbon_chemistry, pH)

function CarbonDioxideConcentration(FT = Float64;
                                    carbon_chemistry = CarbonChemistry(FT),
                                    DIC = :DIC,
                                    Alk = :Alk,
                                    warm_start = false,
                                    grid = nothing)

    pH = warm_start ? warm_start_pH(grid) : nothing

    return CarbonDioxideConcentration{DIC, Alk}(carbon_chemistry, pH)
end

warm_start_pH(::Nothing) = throw(ArgumentError("`warm_start = true` needs a `grid` to store the surface pH on"))
warm_start_pH(grid) = Field{Center, Center, Nothing}(grid)

Adapt.adapt_structure(to, cc::CarbonDioxideConcentration{DIC, Alk}) where {DIC, Alk} =
    CarbonDioxideConcentration{DIC, Alk}(adapt(to, cc.carbon_chemistry), adapt(to, cc.pH))

summary(::CarbonDioxideConcentration{DIC, Alk, CC}) where {DIC, Alk, CC} =
    "`CarbonChemistry` derived aqueous carbon dioxide concentration ([CO₂(aq)], mmol/m³) {$DIC, $Alk, $(nameof(CC))}"

show(io::IO, ccc::CarbonDioxideConcentration{DIC, Alk}) where {DIC, Alk} =
    println(io, summary(ccc), "\n",
            "    Solves the $(nameof(typeof(ccc.carbon_chemistry))) based on $DIC and $Alk",
            isnothing(ccc.pH) ? "" : ", warm started from the stored surface pH")

@inline function surface_value(cc::CarbonDioxideConcentration{DIC_name, Alk_name}, i, j, grid, clock, model_fields) where {DIC_name, Alk_name}
    DIC = @inbounds model_fields[DIC_name][i, j, grid.Nz] # this is a compile time inference so is fine on GPU
    Alk = @inbounds model_fields[Alk_name][i, j, grid.Nz]

    T = @inbounds model_fields.T[i, j, grid.Nz]
    S = @inbounds model_fields.S[i, j, grid.Nz]

    silicate  = silicate_concentration(grid, i, j, grid.Nz, model_fields)
    phosphate = phosphate_concentration(grid, i, j, grid.Nz, model_fields)

    return surface_CO₂(cc.pH, cc.carbon_chemistry, i, j; DIC, Alk, T, S, silicate, phosphate)
end

@inline surface_CO₂(::Nothing, carbon_chemistry, i, j; kwargs...) = carbon_chemistry(; kwargs..., output = Val(:CO₂))

@inline function surface_CO₂(pH, carbon_chemistry, i, j; DIC, kwargs...)
    pH⁻ = @inbounds pH[i, j, 1]

    # nothing stored yet (the field starts at zero) or a bad value falls back to the default guess
    initial_pH_guess = ifelse((pH⁻ > 0) & (pH⁻ < 14), pH⁻, convert(typeof(DIC), 8))

    CO₂, pHⁿ = carbon_chemistry(; DIC, kwargs..., initial_pH_guess, output = Val((:CO₂, :pHᶠ)))

    # each thread only writes its own column so this is safe on the GPU
    @inbounds pH[i, j, 1] = pHⁿ

    return CO₂
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
mmol / m³ per μatm using the TEOS-10 density, which pairs with a
[`CarbonDioxideConcentration`](@ref) water side. Note that this is the solubility of a
*dry air mole fraction* (``f_f = K_0 (1 - p_{H_2O}) \\gamma``), and so already carries the
water vapour and non-ideality corrections; it is not [`K0`](@ref OceanBioME.Models.CarbonChemistryModel.K0).

[`CarbonDioxideGasExchangeBoundaryCondition`](@ref) instead supplies a `solubility` built
from its `carbon_chemistry`'s own `density_function`, so that the air and water sides always
use the same density even if a non-default one is given.

A `solubility` of `nothing` leaves the air concentration as a mole fraction in ppmv.

This drops into the `air_concentration` of a
[`CarbonDioxideGasExchangeBoundaryCondition`](@ref) in place of a bare number, and with the
default pressure of 1 atm it leaves the mole fraction unscaled by pressure.

Note that this is the *only* atmospheric pressure on the carbon dioxide exchange path. The
water-side [`CarbonDioxideConcentration`](@ref) carries no atmospheric pressure (it is evaluated
at the surface, i.e. zero `water_pressure`, and its `[CO₂(aq)]` output does not depend on the
`atmospheric_pressure` of the `CarbonChemistry` call, which only enters its `Val(:pCO₂)` output).
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
                                  solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(FF{FT}(), teos10_polynomial_approximation))

Returns the air-side carbon dioxide concentration ``x_{CO_2} p_{atm} f_f \\rho / 10^3`` (in
mmol / m³), or ``x_{CO_2} p_{atm}`` (in ppmv ≡ μatm) when `solubility = nothing`.

Keyword Arguments
=================

- `mole_fraction`: the dry-air mole fraction of carbon dioxide (ppmv), which may be a
  number, a function of the form `(x, y, t)`, or a `Field`
- `atmospheric_pressure`: the total atmospheric pressure (atm), which may be a number, a
  function of the form `(x, y, t)`, or a `Field`; the default of 1 leaves the mole fraction
  unchanged. Note that the units are atmospheres, not pascals, and that the value is not
  checked
- `solubility`: a function of `(T, S)` returning the conversion from a partial pressure in
  μatm to a concentration in mmol / m³. Defaults to the Weiss and Price (1980) [`FF`](@ref)
  fit with the TEOS-10 density; pass `nothing` to leave the air concentration as a mole
  fraction in ppmv instead. `CarbonDioxideGasExchangeBoundaryCondition` supplies its own
  `carbon_chemistry`'s `density_function` here, so the air and water sides stay matched even
  when a non-default density is used

See also [`CarbonDioxideConcentration`](@ref) and
[`CarbonDioxideGasExchangeBoundaryCondition`](@ref).
"""
function CarbonDioxideAirConcentration(FT = Float64;
                                       mole_fraction = 413,      # ppmv
                                       atmospheric_pressure = 1, # atm
                                       solubility = MolPerKgPerAtmToMMolPerCubicMPerMicroAtm(FF{FT}(), teos10_polynomial_approximation))

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

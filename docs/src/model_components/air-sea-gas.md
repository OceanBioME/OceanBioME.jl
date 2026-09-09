# [Air-sea gas exchange](@id air-sea-gas)

Air-sea gas transfer is typically parameterised as a function of temperature (``T``) and wind speed (``u_{10}``), and the concentration of the gas in the air (``C_a``) and in the surface water (``C_w``) in the form:
```math
F = k(u_{10}, T)(C_w - C_a),
```
where `k` is the gas transfer velocity.

Our implementation is intended to be generic for any gas, so you can specify `air_concentration`, `water_concentration`, `transfer_velocity`, and `wind_speed` as any function in `GasExchange`, but we also provide constructors and default values for carbon dioxide and oxygen.

To setup carbon dioxide and/or oxygen boundary conditions you simply build the condition and then specify it in the model:
```@example gasexchange
using OceanBioME
CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition()
O₂_flux  = OxygenGasExchangeBoundaryCondition()
using Oceananigans

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

model = NonhydrostaticModel(grid;
                            biogeochemistry = LOBSTER(grid; inorganic_carbon = CarbonateSystem(), oxygen = Oxygen()),
                            boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),
                                                    O₂ = FieldBoundaryConditions(top =  O₂_flux)),
                            tracers = (:T, :S))
```

!!! note

    All gas exchange models require temperature (`T`) to be present in the model, and carbon dioxide requires sailinity (`S`), total inorganic carbon (`DIC`), and alkalinity (`Alk`), and optionally can take silicate and phosphate where there names are specified in the keyword argument `silicate_and_phosphate_names`

## Model equations

### Gas transfer velocity

The default gas transfer velocity (`ScaledTransferVelocity`) returns a velocity in the form:
```math
k(u_{10}, T) = cu_{10}^2\left(\frac{Sc(T)}{660}\right)^{-1/2},
```
where ``c`` is a coefficient (`coeff`) which typically is wind product specific with default value ``0.266`` cm/hour from [Ho2006](@citet), and ``Sc`` is gas specific the temperature dependent Schmidt number (the dimensionless ratio of momentum and mass diffusivity) specified as `schmidt_number` which can be any function of temperature. The default parameterisations is the 4th order polynomial formulation of [Wanninkhof2014](@citet).

Currently, the parameters for CO₂ and oxygen are included, but it would be very straightforward to add the parameters given in the original publication for other gases (e.g. inert tracers of other nutrients such as N₂).

### Carbon dioxide concentration

For most gasses the water concentration `C_w` is simply taken directly from the biogeochemical model or another tracer (in which case `water_concentration` should be set to `TracerConcentration(:tracer_name)`), but for carbon dioxide it must be derived from the dissolved inorganic carbon (`DIC`) and `Alk`alinity by a `CarbonChemistry` model (please see the docs for [CarbonChemistry](@ref carbon-chemistry)).

The water concentration is the aqueous carbon dioxide concentration in mmol / m³,
```math
C_w = [CO_2(aq)] = DIC\frac{[H^+]^2}{[H^+]^2 + K_1[H^+] + K_1K_2},
```
(`CarbonDioxideConcentration`), and the air concentration is the dry air mole fraction converted onto the same basis by Dalton's law and a solubility,
```math
C_a = x(CO_2)p_{atm}f_f(T, S)\frac{\rho}{10^3}.
```
The solubility ``f_f`` is the [Weiss1980](@citet) parameterisation (`FF`), which is the solubility ``K_0`` corrected for the water vapour pressure of saturated air and for the non-ideality of the gas phase,
```math
f_f = K_0(1 - p_{H_2O})\gamma,
```
and is therefore the right quantity to multiply a *dry air* mole fraction by. The transfer velocity is a bare piston velocity (its `solubility` is `UnitSolubility`), and the atmospheric pressure enters the flux exactly once, on the air side.

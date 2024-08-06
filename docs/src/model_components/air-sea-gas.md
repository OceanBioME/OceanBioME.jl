# [Air-sea gas exchange](@id air-sea-gas)

Air-sea gas transfer is typically parameterised as a function of temperature (``T``) and wind speed (``u_{10}``), and the concentration of the gas in the air (``C_a``) and in the surface water (``C_w``) in the form:
```math
F = k(u_{10}, T)(C_w - C_a),
```
where `k` is the gas transfer velocity.

Our implementation is intended to be generic for any gas so you can specify `air_concentration`, `water_concentration`, `transfer_velocity`, and `wind_speed` as any function in `GasExchange`, but we also provide constructors and default values for carbon dioxide and oxygen. 

To setup carbon dioxide and/or oxygen boundary conditions you simply build the condition and then specify it in the model:
```@example gasexchange
using OceanBioME
CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition()
O₂_flux  = OxygenGasExchangeBoundaryCondition()
using Oceananigans

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

model = NonhydrostaticModel(; grid,
                              biogeochemistry = LOBSTER(; grid, carbonates = true, oxygen = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), 
                                                      O₂ = FieldBoundaryConditions(top =  O₂_flux)),
                              tracers = (:T, :S))
```

!!! compat Field dependencies

    All gas exchange models require temperature (`T`) to be present in the model, and carbon dioxide requires sailinity (`S`), total inorganic carbon (`DIC`), and alkalinity (`Alk`), and optionally can take silicate and phosphate where there names are specified in the keyword argument `silicate_and_phosphate_names`

## Model equations

### Gas transfer velocity

The default gas transfer velocity (`ScaledTransferVelocity`) returns a velocity in the form:
```math
k(u_{10}, T) = cu_{10}^2\left(\frac{Sc(T)}{660}\right)^{-1/2},
```
where ``c`` is a coefficient (`coeff`) which typically is wind product specific with default value ``0.266`` cm/hour from [Ho2006](@citet), and ``Sc`` is gas specific the temperature dependent Schmidt number (the dimensionless ratio of momentum and mass diffusivity) specified as `schmidt_number` which can be any function of temperature. The default parameterisations is the 4th order polynomial formulation of [Wanninkhof2014](@citet).

Currently, the parameters for CO₂ and oxygen are included, but it would be very straightforward to add the parameters given in the original publication for other gases (e.g. inert tracers of other nutrients such as N₂).

### Carbon dioxide partial pressure

For most gasses the water concentration `C_w` is simply taken directly from the biogeochemical model or another tracer (in which case `water_concentration` should be set to `TracerConcentration(:tracer_name)`), but for carbon dioxide the fugacity (``fCO_2``) must be derived from the dissolved inorganic carbon (`DIC`) and `Alk`alinity by a `CarbonChemistry` model (please see the docs for [CarbonChemistry](@ref carbon-chemistry)), and used to calculate the partial pressure (``pCO_2``).

The default parameterisation for the partial pressure (`CarbonDioxideConcentration`) is given by [dickson2007](@citet) and defines the partial pressure to be the mole fraction ``x(CO_2)`` multiplied by the pressure, ``P``, related to the fugacity by:
```math
fCO_2 = x(CO_2)P\exp\left(\frac{1}{RT}\int_0^P\left(V(CO_2)-\frac{RT}{P'}\right)dP'\right).
```
The volume (``V``) is related to the gas pressure by the virial expression:
```math
\frac{PV(CO_2)}{RT}\approx1+\frac{B(x, T)}{V(CO_2)}+\mathcal{O}(V(CO_2)^{-2}),
```
and the first virial coefficient ``B`` for carbon dioxide in air can be approximated as:
```math
B_{CO_2-\text{air}} \approx B_{CO_2}(T) + 2x(CO_2)\delta_{CO_2-\text{air}}(T),
```
where ``\delta`` is the cross virial coefficient.

``B_{CO_2}`` and ``\delta_{CO_2-\text{air}}`` are parameterised by [Weiss1974](@citet) and reccomended in [dickson2007](@citet) as fourth and first order polynomials respectively.

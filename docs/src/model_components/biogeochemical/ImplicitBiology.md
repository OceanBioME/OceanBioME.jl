# [ImplicitBiology](@id ImplicitBiology)

`ImplicitBiology` computes a community productivity limited by nutrient availability and light, but does not explicitly track any planktonic biomass. This makes it suitable for large-scale or long-timescale simulations where resolving full plankton dynamics is unnecessary or too expensive. The formulation follows the Darwin/MITGCM implicit biology approach of [dutkiewicz_sokolov_scott_stone_2005](@citet).

The default tracer set is `(:PO₄, :DOP, :POP, :DIC, :Alk)`, though this varies with the choice of limiting nutrients and detritus model.

```@raw html
<!-- placeholder for diagram -->
```

## Model equations

Community productivity is parameterised as:

```math
P = \alpha_\text{max} L_N L_{PAR} L_T,
```

where ``\alpha_\text{max}`` is the maximum community productivity, ``L_T`` is an optional temperature dependence (one by default), and the nutrient and light limitations are:

```math
L_{PAR} = \frac{PAR}{PAR + k_{PAR}},
```

```math
L_N = \min\left(\frac{PO_4}{PO_4 + k_{PO_4}},\ \frac{NO_3}{NO_3 + k_{NO_3}},\ \frac{Fe}{Fe + k_{Fe}}\right).
```

The productivity drives nutrient uptake and waste production. The waste is partitioned between a dissolved and a particulate pool:

```math
\frac{\partial DOP}{\partial t} = f_d P - \mu_{DOP} DOP,
```

```math
\frac{\partial POP}{\partial t} = (1 - f_d) P - \mu_{POP} POP,
```

where ``f_d`` is the dissolved fraction of waste, and ``\mu_{DOP}`` and ``\mu_{POP}`` are remineralisation rates. Nutrients are replenished by remineralisation of both pools:

```math
\frac{\partial PO_4}{\partial t} = -P + \mu_{DOP} DOP + \mu_{POP} f_{r} POP,
```

where ``f_r`` is the dissolved fraction of particulate remineralisation (defaulting to 1, i.e. remineralisation returns fully to the dissolved inorganic pool). Nitrogen and iron follow the same structure scaled by their stoichiometric ratios ``r_N = N:P`` and ``r_{Fe} = Fe:P``.

When `inorganic_carbon = CarbonateSystem()` (the default), dissolved inorganic carbon and alkalinity evolve as described in the [LOBSTER carbonate chemistry section](@ref LOBSTER).

### Parameter variable names

| Symbol                    | Variable name                       | Units           |
|---------------------------|-------------------------------------|-----------------|
| ``\alpha_\text{max}``     | `maximum_community_productivity`    | mmol P / m³ / s |
| ``k_{PAR}``               | `light_half_saturation`             | W / m²          |
| ``f_d``                   | `dissolved_fraction_of_waste`       | -               |
| ``r_C``                   | `carbon_ratio`                      | mol C / mol P   |
| ``r_N``                   | `nitrogen_ratio`                    | mol N / mol P   |
| ``r_{Fe}``                | `iron_ratio`                        | mol Fe / mol P  |
| ``\rho_{CaCO_3}``         | `rain_ratio`                        | mol CaCO₃ / mol C |
| ``k_{PO_4}``              | `nutrient_half_saturations.phosphate` | mmol P / m³   |
| ``k_{NO_3}``              | `nutrient_half_saturations.nitrate`   | mmol N / m³   |
| ``k_{Fe}``                | `nutrient_half_saturations.iron`      | mmol Fe / m³  |
| ``L_T``                   | `temperature_dependence`              | -             |

Default detritus parameter values are given in [Parameters](@ref parameters).

## Setup

The default configuration limits productivity by nitrate, phosphate, and iron:

```julia
using OceanBioME, Oceananigans

grid = RectilinearGrid(size = 10, extent = 200, topology = (Flat, Flat, Bounded))

model = NonhydrostaticModel(grid; biogeochemistry = ImplicitBiology(grid))
```

The set of limiting nutrients can be changed with `limiting_nutrients`:

```julia
# Limit by nitrate only
biogeochemistry = ImplicitBiology(grid; limiting_nutrients = (:nitrate,))

# Limit by nitrate and iron (omits phosphate)
biogeochemistry = ImplicitBiology(grid; limiting_nutrients = (:nitrate, :iron))
```

Productivity is independent of temperature by default. An exponential [Eppley (1972)](@cite eppley_1972) response, ``L_T = Q_{10}^{(T - T_r)/10}`` with ``Q_{10} = 1.88`` (Eppley's ``1.066^T``) and ``T_r = 0``°C by default, is added with `EppleyCurve`, which then requires a temperature tracer `T`. ``\alpha_\text{max}`` is then the productivity at ``T_r``:

```julia
biogeochemistry = ImplicitBiology(grid; temperature_dependence = EppleyCurve())

# a Q₁₀ of 2 referenced to 20°C
biogeochemistry = ImplicitBiology(grid; temperature_dependence = EppleyCurve(temperature_coefficient = 2, reference_temperature = 20))
```

Any object callable on the temperature (°C) can be supplied instead of `EppleyCurve` to define a different response.

Carbonate chemistry and oxygen can be added as with any other NPD model:

```julia
biogeochemistry = ImplicitBiology(grid; oxygen = Oxygen())
```

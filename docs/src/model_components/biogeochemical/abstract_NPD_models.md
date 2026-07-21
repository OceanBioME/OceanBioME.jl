# [Nutrients-Plankton-Detritus framework](@id npd_framework)

The Nutrient-Plankton-Detritus framework allows biogeochemical models to be built from shared components. Individual biogeochemical models like [LOBSTER](@ref LOBSTER) or [NPZD](@ref NPZD) can be thought of as 'recipes' that combine these ingredients in different ways. 

The framework includes five component types: Nutrients, Plankton (including phytoplankton and zooplankton), Detritus, Inorganic Carbon, and Oxygen. This framework is extremely flexible. Models built using this framework don't need to use all of these component types, but they can also include mutiple components. For example, multiple nutrient and phytoplankton types can be used in the biogeochemical model and it is even possible not to include explicit phytoplankton (as in [ImplicitBiology](@ref ImplicitBiology)). This framework allows biogeochemical models to be highly composable, and greatly aids development of new biogeochemical models. For a step-by-step guide to creating a new biogeochemical model using the Nutrient-Plankton-Detritus framework, see [Implementing new models](@ref model_implementation).

The `NutrientsPlanktonDetritus` framework organises the biogeochemistry into five component types. Biogeochemical models must specify which of the existing types are used for each component (or use the default of `nothing` in the case of Inorganic Carbon and Oxygen).

| Component | Keyword | Available types |
|------|---------|-----------------|
| Nutrients | `nutrients` | `Nutrients` (constructed with `NitrateAmmonia`, `N`, `PO₄`, `Fe`, or `Si` sub-components) |
| Plankton | `plankton` | `Abiotic`, `ImplicitProductivity`, `PhytoZoo` |
| Detritus | `detritus` | `InstantRemineralisationDetritus`, `Detritus`, `DissolvedParticulate`, `CarbonNitrogenDissolvedParticulate` |
| Inorganic Carbon | `inorganic_carbon` | `nothing` (default), `CarbonateSystem` |
| Oxygen | `oxygen` | `nothing` (default), `Oxygen` |

## Component structure

Components fall into two categories: **configuration objects** that hold scalar parameters and specify which tracers to include, and **spatial components** that may also hold grid-aligned fields (for example, sinking velocities). The `grid` argument to `NutrientsPlanktonDetritus(grid; ...)` is always required to build the full model; individual components only need `grid` when they configure such fields.

### Nutrients

[`Nutrients`](@ref) groups the inorganic nutrient pools that can limit plankton growth. It has four fixed slots — nitrogen, phosphate, iron, and silicate — each of which may be `nothing` (that nutrient is not tracked) or one of the types below. `Nutrients` stores only scalar parameters and does not require `grid`; the actual tracer fields are created when the model is assembled.

When using positional arguments, the four slots are passed in order: `Nutrients(nitrogen, phosphate, iron, silicate)`. The keyword `nothing` can be used to omit any of the nutrients. For example, `Nutrients(NitrateAmmonia(), PO₄, Fe, nothing)` tracks nitrate, ammonia, phosphate, and iron, but not silicate.

- **`N`** — a single combined nitrogen tracer (`N`).
- **`NitrateAmmonia`** — separate nitrate and ammonia tracers (`NO₃`, `NH₄`) with nitrification between them.
- **`PO₄`** — a phosphate tracer.
- **`Fe`** — an iron tracer.
- **`Si`** — a silicate tracer.

### Plankton

Phytoplankton and zooplankton are grouped together in the `Plankton` component. The available component types for building biogeochemical models using the NPD framework are:

- **`Abiotic`** (default) — phytoplankton and zooplankton are not explicitly tracked. Adds no tracers and produces no biological source or sink terms. Useful for purely abiotic experiments (for example, carbonate chemistry without biology).
- **`ImplicitProductivity`** — nutrient- and light-limited community productivity with no explicit plankton biomass (used by [`ImplicitBiology`](@ref ImplicitBiology)).
- **`PhytoZoo`** — explicit phytoplankton (`P`) and zooplankton (`Z`) with growth, grazing, and loss terms (used by [`LOBSTER`](@ref LOBSTER) and [`NPZD`](@ref NPZD)). Call `PhytoZoo(grid; ...)` to configure sinking-velocity fields on the grid; call `PhytoZoo(; ...)` without `grid` when sinking is not needed (for example, in a box model).

### Detritus

There are several options for handling detritus in the NPD constructor. Detritus can be treated implicitly by assuming that organic waste is instantly remineralised and made available in the dissolved inorganic nutrient pool. Detritus in dissolved and particulate forms can also be tracked explicitly with varying degrees of complexity.

- **`InstantRemineralisationDetritus`** (default) — detritus is not explicitly tracked as a tracer. Instead, the model returns plankton waste straight into the inorganic nutrient pool. Useful for closed-budget box models.
- **`Detritus`** — a single sinking detritus pool (`D`); requires `grid` as an argument to configure sinking speeds.
- **`DissolvedParticulate`** — configurable dissolved and particulate organic pools (for example, `DOM`, `sPOM`, `bPOM`); requires `grid` to configure particulate sinking speeds.
- **`CarbonNitrogenDissolvedParticulate`** — variable-Redfield carbon and nitrogen tracked separately in dissolved, small-particulate, and large-particulate pools (`DON`, `DOC`, `sPON`, `sPOC`, `bPON`, `bPOC`); requires `grid`.

### Inorganic carbon

The `Inorganic Carbon` component can be used to enable carbon chemistry. This component can be used with or without explicit biology.

- **`nothing`** (default) — inorganic carbon is not tracked.
- **`CarbonateSystem`** — dissolved inorganic carbon and alkalinity (`DIC`, `Alk`), driven by primary production, remineralisation, and implicit calcite production/dissolution.

### Oxygen

The `Oxygen` component can be used to track dissolved oxygen.

- **`nothing`** (default) — oxygen is not tracked.
- **`Oxygen`** — an oxygen tracer (`O₂`) produced by photosynthesis and consumed by remineralisation and nitrification. Oxygen is one-way coupled, meaning that the current implementation of biology does not depend on oxygen.


Users can over-ride the default parameters for any component by optionally passing in parameter values. Not all parameters need to be specified.  Any parameters that are not specified will be set to the default value.


## Preset constructors

For common configurations, convenience constructors assemble the components using the 'recipe' for that biogeochemical model with default parameters. You can build each biogeochemical model with default parameters simply by passing the model grid as an argument:

- `LOBSTER(grid)` — medium-complexity model with phytoplankton, zooplankton, nitrate, ammonia, and dissolved/particulate detritus. See the [LOBSTER](@ref LOBSTER) page.
- `NPZD(grid)` — simple four-compartment nutrient–phytoplankton–zooplankton–detritus model from [Kuhn2015](@citet). See the [NPZD](@ref NPZD) page.
- `ImplicitBiology(grid)` — nutrient-limited community productivity with no explicit plankton biomass. See the [ImplicitBiology](@ref ImplicitBiology) page.

Any of the default parameter values can be over-ridden by passing in parameter values to the NPD components. The following example builds a NPZD model with user-specified phytoplankton and zooplankton growth rates:

```julia
using OceanBioME, Oceananigans

grid = RectilinearGrid(size = 10, extent = 200, topology = (Flat, Flat, Bounded))

model = NonhydrostaticModel(grid;
                            biogeochemistry = NPZD(grid;
                                                   plankton = PhytoZoo(grid;
                                                                      phytoplankton_maximum_growth_rate = 1 / day,
                                                                      maximum_grazing_rate = 2 / day)))
```

All preset constructors accept the same optional keyword arguments as `NutrientsPlanktonDetritus` itself (e.g.`inorganic_carbon`, `oxygen`, `light_attenuation`, `sediment`, and `scale_negatives`) so you can extend any of the constructors without building a model from scratch. For example, the following code adds inorganic carbon and oxygen to the LOBSTER model:

```julia
using OceanBioME, Oceananigans

grid = RectilinearGrid(size = 10, extent = 200, topology = (Flat, Flat, Bounded))

model = NonhydrostaticModel(grid;
                            biogeochemistry = LOBSTER(grid;
                                                      inorganic_carbon = CarbonateSystem(),
                                                      oxygen = Oxygen()))
```

## Custom configurations

You can create a fully custom biogeochemical model by calling `NutrientsPlanktonDetritus` directly and specifying each component. For example, to build a model with Nitrate, Ammonia, Phosphate, and Iron, one phytoplankton and one zooplankton type, one dissolved detritus pool, and an active carbonate system: 

```julia
biogeochemistry = NutrientsPlanktonDetritus(grid;
                                            nutrients  = Nutrients(NitrateAmmonia(), PO₄, Fe, nothing),
                                            plankton   = PhytoZoo(grid),
                                            detritus   = DissolvedParticulate(grid),
                                            inorganic_carbon = CarbonateSystem())
```

The framework dispatches through the components to assemble the required tracers and auxiliary fields automatically, so all standard Oceananigans boundary conditions and output writers work as normal.


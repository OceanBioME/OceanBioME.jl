# [Western Tropical South Pacific (WTSP) diazotroph model](@id WTSP)

The `WTSP` model is a six-tracer model of an oligotrophic, iron limited surface ocean in which
dinitrogen fixation, rather than the supply of fixed nitrogen, sets new production. It was built for
the Western Tropical South Pacific, where `Trichodesmium`-like diazotrophs bloom on iron delivered by
shallow hydrothermal and volcanic sources, but nothing in it is specific to that region.

It is a preset of the [Nutrients-Plankton-Detritus framework](@ref npd_framework), assembled from
[`Nutrients`](@ref) (nitrogen, phosphate, and iron), a [`PhytoDiazotroph`](@ref) plankton, and a
single-class [`Detritus`](@ref) pool:

```julia
using OceanBioME, Oceananigans

grid = RectilinearGrid(size = 10, extent = 200, topology = (Flat, Flat, Bounded))

model = NonhydrostaticModel(grid; biogeochemistry = WTSP(grid))
```

which evolves the tracers `N`, `PO₄`, `Fe`, `P`, `Diaz`, and `D`. As with the other presets, any
component can be replaced or reparameterised, and [`CarbonateSystem`](@ref) and [`Oxygen`](@ref) can
be added.

The distinguishing feature is the second photoautotroph. Picoplankton (`P`) grow on dissolved
inorganic nitrogen, phosphate, iron, and light in the usual way. Diazotrophs (`Diaz`) need no fixed
nitrogen at all, but their nitrogenase makes them iron hungry — their iron quota is about eight times
the picoplankton's, and their iron half-saturation about twenty-five times larger — so iron
availability, not nitrogen, decides where they can live. With the default parameters they need
dissolved iron above about 2 nmol / m³ to outgrow their mortality, well above the ~0.3 nmol / m³ of
the surrounding gyre.

Fixed nitrogen reaches the rest of the ecosystem along two paths: a fraction ``\gamma`` of gross
fixation is released directly as dissolved inorganic nitrogen, and the rest enters biomass and is
returned when diazotrophs die and their detritus is remineralised.

## Model equations

Writing ``\mathrm{DIN}`` for the dissolved inorganic nitrogen (`N`, or `NO₃` + `NH₄`), each group has
a light and Liebig-nutrient limited specific growth rate

```math
\mu_P = \mu_P^\text{max}\left(1 - e^{-PAR/k_{PAR}^P}\right)\min\left(\frac{\mathrm{DIN}}{\mathrm{DIN} + k_N},\ \frac{PO_4}{PO_4 + k_{PO_4}^P},\ \frac{Fe}{Fe + k_{Fe}^P}\right),
```

```math
\mu_D = \mu_D^\text{max}\left(1 - e^{-PAR/k_{PAR}^D}\right)\min\left(\frac{PO_4}{PO_4 + k_{PO_4}^D},\ \frac{Fe}{Fe + k_{Fe}^D}\right)\frac{k_{\mathrm{DIN}}}{k_{\mathrm{DIN}} + \mathrm{DIN}},
```

where the diazotroph has no nitrogen term but is inhibited by dissolved inorganic nitrogen, which is
energetically cheaper for it to use than fixing N₂. Gross nitrogen fixation is then
``F = \mu_D\,\mathrm{Diaz}``, and losses are ``\ell_X = m_X X + m_X' X^2`` (the density dependent part
is off by default, and closes the model where grazers, which are not resolved, would otherwise be
needed).

```math
\frac{\partial P}{\partial t} = \mu_P P - \ell_P,
```

```math
\frac{\partial \mathrm{Diaz}}{\partial t} = (1 - \gamma) F - \ell_{\mathrm{Diaz}},
```

```math
\frac{\partial D}{\partial t} = \ell_P + \ell_{\mathrm{Diaz}} - r D,
```

```math
\frac{\partial \mathrm{DIN}}{\partial t} = \gamma F + r D - \mu_P P,
```

```math
\frac{\partial PO_4}{\partial t} = R_P\left(r D - \mu_P P - (1 - \gamma) F\right),
```

```math
\frac{\partial Fe}{\partial t} = R_{Fe}\left(r D - \mu_P P - (1 - \gamma) F\right) + \Delta_{Fe}\left(\ell_{\mathrm{Diaz}} - (1 - \gamma)F\right),
```

where ``R_P`` and ``R_{Fe}`` are the phosphorus and iron content of all organic matter other than the
diazotroph's extra iron, and ``\Delta_{Fe} = R_{Fe}^{\mathrm{Diaz}} - R_{Fe}`` is that extra quota.
Detritus additionally sinks, and phytoplankton may be given a sinking (or, for `Trichodesmium`,
rising) speed.

Only the nitrogen released by fixation and the nitrogen fixed into diazotroph biomass enter the model
from outside; the ``\gamma F`` term carries no carbon, phosphorus, or iron, because it comes straight
from N₂.

## Model conservation

Nitrogen fixation is the only source of any element, so

```math
\frac{\partial}{\partial t}\left(\mathrm{DIN} + P + \mathrm{Diaz} + D\right) = F,
```

and phosphorus, iron, carbon, and — when [`Oxygen`](@ref) is included — oxygen are conserved exactly
(excluding sinking and external sources). Setting `Diaz` to zero removes the source and nitrogen is
conserved too.

Iron is the one element whose inventory is not a single ratio times the tracers, since the diazotroph
carries a larger quota than the detritus its dead cells become. The excess ``\Delta_{Fe}`` is
therefore returned straight to the dissolved pool as diazotrophs die, rather than being buried in a
detritus pool which cannot represent two quotas at once. [`conserved_tracers`](@ref) reports the
weight each tracer carries, which for iron is

```math
Fe + R_{Fe}\left(P + D\right) + R_{Fe}^{\mathrm{Diaz}}\,\mathrm{Diaz}.
```

Two framework-level terms make the carbon and oxygen budgets right in the presence of fixation:
alkalinity gains the phosphate which diazotrophs take up alongside nitrogen they did not take from
the dissolved inorganic nitrogen pool, and oxygen is not credited the nitrate-reduction oxygen for
nitrogen which was fixed from N₂ rather than reduced from nitrate. Both are zero for every plankton
which does not fix nitrogen.

Iron scavenging, atmospheric and hydrothermal iron supply, lateral nutrient supply, and any other
external source are deliberately *not* part of the biogeochemistry. They are properties of a
particular setup rather than of the ecosystem, and are best added as Oceananigans forcings or
boundary conditions, which also keeps the elemental budgets above exactly diagnosable.

## Parameter variable names

Rates are stored per second; the values below are given per day where that is more natural.

| Symbol | Variable name | Value | Units |
|--------|---------------|-------|-------|
| $\mu_P^\text{max}$ | `picoplankton_maximum_growth_rate` | 0.8 | 1 / day |
| $\mu_D^\text{max}`` | `diazotroph_maximum_growth_rate` | 0.2 / 0.7 | 1 / day |
| $k_N$ | `picoplankton_nutrient_half_saturations.nitrate` | 0.1 | mmol N / m³ |
| $k_{PO_4}^P$ | `picoplankton_nutrient_half_saturations.phosphate` | 0.03 | mmol P / m³ |
| $k_{Fe}^P$ | `picoplankton_nutrient_half_saturations.iron` | 8 × 10⁻⁵ | mmol Fe / m³ |
| $k_{PO_4}^D$ | `diazotroph_nutrient_half_saturations.phosphate` | 0.05 | mmol P / m³ |
| $k_{Fe}^D$ | `diazotroph_nutrient_half_saturations.iron` | 2 × 10⁻³ | mmol Fe / m³ |
| $k_{PAR}^P$, ``k_{PAR}^D`` | `picoplankton_light_half_saturation`, `diazotroph_light_half_saturation` | 2.17 | W / m² |
| $k_{\mathrm{DIN}}$ | `nitrogen_fixation_inhibition_half_saturation` | 0.3 | mmol N / m³ |
| $\gamma$ | `diazotroph_nitrogen_release_fraction` | 0.3 | - |
| $m_P$ | `picoplankton_mortality_rate` | 0.06 | 1 / day |
| $m_{\mathrm{Diaz}}$ | `diazotroph_mortality_rate` | 0.1 | 1 / day |
| $m_P'$, $m_{\mathrm{Diaz}}'$ | `picoplankton_quadratic_mortality_rate`, `diazotroph_quadratic_mortality_rate` | 0 | 1 / day / mmol N / m³ |
| $r$ | `remineralisation_rate` ([`Detritus`](@ref)) | 0.4 | 1 / day |
| $w_D$ | `sinking_speed` ([`Detritus`](@ref)) | 1.5 | m / day |
| $R_C$ | `carbon_ratio` | 6.625 | mol C / mol N |
| $R_P$ | `phosphate_ratio` | 1/16 | mol P / mol N |
| $R_{Fe}$ | `picoplankton_iron_ratio` | 3.3125 × 10⁻⁵ | mol Fe / mol N |
| $R_{Fe}^{\mathrm{Diaz}}$ | `diazotroph_iron_ratio` | 2.65 × 10⁻⁴ | mol Fe / mol N |

The two iron quotas correspond to Fe:C of 5 and 40 μmol Fe / mol C, and the chlorophyll ratio of
1.137 g Chl / mol N to a C:Chl of 70 g C / g Chl.

Note that $\mu_D^\text{max}$ is *gross* of the released fraction, so diazotroph biomass grows at
most at $(1 - \gamma)\mu_D^\text{max} = 0.2$/day.

## Adding external supply

An oligotrophic box needs some supply to reach a steady state. Here iron is restored towards a target
concentration, standing in for a shallow hydrothermal or volcanic source, while everything else stays
in a closed budget:

```julia
using OceanBioME, Oceananigans, Oceananigans.Units
using Oceananigans.Fields: ConstantField

grid = BoxModelGrid()
clock = Clock(time = zero(grid))

biogeochemistry = WTSP(grid;
                       light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(70)))

restore(name, rate, target) = (clock, fields) -> rate * (target - @inbounds fields[name][1, 1, 1])

dFe = 3e-3 # mmol Fe / m³

model = BoxModel(; biogeochemistry, grid, clock,
                   forcing = (; Fe = restore(:Fe, 1/30days, dFe)))

set!(model, N = 0.05, PO₄ = 0.15, Fe = dFe, P = 0.05, Diaz = 1e-4, D = 0.01)

simulation = Simulation(model; Δt = 10minutes, stop_time = 3 * 365days)

run!(simulation)
```

Sweeping `dFe` shows the iron threshold: at 0.3 and 1 nmol / m³ the diazotrophs seeded at
`Diaz = 1e-4` die out and phosphate stays near its initial 0.15 mmol / m³, while at 3 and 6 nmol / m³
they bloom, draw phosphate down to 0.01–0.02 mmol / m³, and raise picoplankton biomass by more than
an order of magnitude. The bloom is not a steady state: it runs on the standing stock of phosphate,
and once that is drawn down the diazotrophs die back again.

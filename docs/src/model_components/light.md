# [Light attenuation models](@id light)

Nearly all BGC models require some model of the attenuation of PAR through the water. Usually this depends on the concentration of chlorophyll in the water (in phytoplankton), and may depend on the concentration of coloured dissolved organic matter or particulates.

We have two models implemented, a two band model by [Karleskind2011](@citet), and a more generic "multi band" model which can have the PAR split into arbitary many wavelength bands, but default to the widely used three band model by [Morel1988](@citet). As the light level is diagnostic of the phytoplankton concentration these models are implemented with the light level as various auxiliary fields which are updated within the biogeochemical model.

Models requiring light attenuation models will set these up automatically, for example [LOBSTER](@ref LOBSTER) sets `light_attenuation_model = TwoBandPhotosyntheticallyActiveRadiation()`. You may choose others. Additionally, you can pass the surface PAR as a function of horizontal position and time. The default for LOBSTER is `(x, y, t) -> 100*max(0.0, cos(t*π/(12hours)))`.

The two band and prescribed attenuation models (below) share a common internal implementation of exponential light attenuation, which also allows them to optionally track ``PAR`` at cell faces as well as cell centres — see [Recording PAR at cell faces](@ref) below.

## The multi band model
The surface intensity is split into multiple bands (usually with equal weight, but users may specify custom weights), and the attenuation of each band (i) is computed from the radiative transfer equation:
```math
\frac{\partial PAR^i}{\partial z} = PAR^i (k^w(i) + \chi(i)Chl^{e(i)}),
```
where ``Chl`` is the concentration of chlorophyll, ``k^w(i)`` is the band specific water attenuation coefficient, ``\chi(i)`` the chlorophyll attenuation coefficient, and ``e(i)`` the chlorophyll exponent.

The water concentration of chlorophyll is returned by a function `chlorophyll` with arguments `biogeochemistry` and `model`. For the `LOBSTER` model this returns a constant ratio of the phytoplankton concentration, but may be different for other models.

## The two band model
Light attenuation is calculated by integrating attenuation (from the surface). The ``PAR`` is considered as two components attenuated at different rates. At depth $z$ the total $PAR$ is given by:

```math
PAR = \frac{PAR_0}{2} \left[\exp\left(k_rz + \chi_r\int_{z=0}^z Chl_r dz\right) + \exp\left(k_bz + \chi_b\int_{z=0}^z Chl_b dz\right)\right],
```

where ``PAR_0`` is the surface value, ``k_r`` and ``k_b`` are the red and blue attenuation coefficients of water, ``\chi_r`` and ``\chi_b`` are the red and blue chlorophyll attenuation coefficients, and ``Chl_r`` and ``Chl_b`` are the red and blue chlorophyll pigment concentrations. Total chlorophyll, ``Chl``, is obtained from `chlorophyll(model.biogeochemistry, model)`. 

 The red and blue pigment concentrations are then found as ``Chl_r = \left(\frac{Chl}{r_\text{pig}}\right)^{e_r}`` and ``Chl_b = \left(\frac{Chl}{r_\text{pig}}\right)^{e_b}``.

### Parameter variable names

| Symbol           | Variable name                     | Units                             |
|------------------|-----------------------------------|-----------------------------------|
| ``k_r``          | `water_red_attenuation`           | 1 / m                             |
| ``k_b``          | `water_blue_attenuation`          | 1 / m                             |
| ``\chi_r``       | `chlorophyll_red_attenuation`     | 1 / m / (mg Chl / m³) ``^ {e_r}`` |
| ``\chi_b``       | `chlorophyll_blue_attenuation`    | 1 / m / (mg Chl / m³) ``^ {e_b}`` |
| ``e_r``          | `chlorophyll_red_exponent`        | -                                 |
| ``e_b``          | `chlorophyll_blue_exponent`       | -                                 |
| ``r_\text{pig}`` | `pigment_ratio`                   | -                                 |

## Recording PAR at cell faces

By default, the light attenuation models above return ``PAR`` sampled at the *centre* of each grid cell, found by integrating the attenuation downward from the surface. `TwoBandPhotosyntheticallyActiveRadiation` and `PrescribedAttenuationPAR` can additionally be given a field on which to record the attenuation at cell *faces* (interfaces), by passing the `interface_field` keyword argument, for example:

```julia
using OceanBioME
using Oceananigans.Fields: ZFaceField

interface_field = ZFaceField(grid)
light_attenuation = TwoBandPhotosyntheticallyActiveRadiation(grid, surface_PAR; interface_field)
```

When `interface_field` is provided, attenuation is instead computed directly between successive interfaces (faces), and the value stored at each cell centre is the exponentially-weighted average of ``PAR`` over that cell, consistent with the tracked face values, rather than a pointwise sample taken at the centre. This is probably a more accurate representation of the light available to a cell than the default pointwise sampling. If ``K`` is the local attenuation coefficient and ``PAR_{k+1}`` is the value at the shallower face of a cell of thickness ``\Delta z``, the cell-averaged value used is

```math
\overline{PAR} = PAR_{k+1} \frac{1 - e^{-K \Delta z}}{K \Delta z}.
```

For `TwoBandPhotosyntheticallyActiveRadiation`, the face values are also exposed as an additional `PAR_interface` auxiliary field. This can be useful on coarse grids, or wherever a flux-consistent estimate of the light available over a whole cell is preferred to a single point sample.

!!! note "Experimental"
    `interface_field` is a new option. For the two band model the red and blue bands' face-to-face attenuations are combined before the cell average above is taken, and this combined-band average has not yet been fully verified — treat it as experimental until confirmed against a reference solution.

## Prescribed attenuation model

`PrescribedAttenuationPAR` is a simpler model that integrates a prescribed attenuation coefficient downward from the surface, rather than computing attenuation from the chlorophyll concentration. At depth ``z`` the PAR is:

```math
PAR(z) = PAR_0 \exp\!\left(-\int_0^z K\, dz'\right),
```

where ``K`` is the attenuation coefficient (units: 1/m). By default ``K`` is a constant, but it may also be a function of position and time.

This model is used by default in `ImplicitBiology` and is useful when no chlorophyll-based feedback on light is desired. It is set up as:

```julia
using OceanBioME

light_attenuation = PrescribedAttenuationPAR(grid, surface_PAR; attenuation = 0.1)
```

where `surface_PAR` may be a constant or a function `f(x, y, t)`, and `attenuation` may be a constant or a function `f(x, y, z, t)`. `PrescribedAttenuationPAR` also accepts the `interface_field` keyword described in [Recording PAR at cell faces](@ref).
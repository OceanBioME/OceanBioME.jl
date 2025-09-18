# [Simple multi-G](@id multi-g)

This model, proposed by [Soetaert2000](@citet), is a "G class" model that evolves carbon and nitrogen in three classes (fast, slow and refectory). The model is also only compatible with the LOBSTER biogeochemical model with carbonate chemistry, oxygen, and variable redfield options on. You also must ensure that the `open_bottom` option is on for particles to leave the bottom of the domain to the sediment model.

It is straightforward to set up:

```jldoctest simplemultig; filter = r".*@ OceanBioME.Models.Sediments.*"
using OceanBioME, Oceananigans, OceanBioME.Sediments

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200))

sediment = SimpleMultiGSediment(grid)

# output
`BiogeochemicalSediment` with `Single-layer multi-G sediment model (Float64)` biogeochemsitry
    Prognostic fields: (:Ns, :Nf, :Nr)
    Tracked fields: (:NO₃, :NH₄, :O₂, :sPOM, :bPOM)
    Coupled fields: (:NO₃, :NH₄, :O₂)
```

You may optionally specify the model parameters. This can then be passed in the setup of a BGC model:

```julia
biogeochemistry = LOBSTER(; grid,
                            carbonate_system = CarbonateSystem(),
                            oxygen = Oxygen(), 
                            detritus = VariableRedfieldDetritus(; open_bottom = true),
                            sediment)
```

### Model equations

This model evolved the carbon and nitrogen components of three liability classes: fast, slow, and refractory. Each component is remineralised with first order decay so evolved like:

```math
\frac{dX_i}{dt} = F_{X_i} - \lambda_iX_i.
```

For the fast and slow classes ``\lambda`` is a positive, non-zero, rate constant, and for the refractory class it is ``0``. ``F_{X_i}`` is the flux, and the flux for each class is simply a constant fraction of the total flux:

```math
F_{X_i} = f_iF_X.
```

The fraction of remineralised sediment (``X_\text{min} = \Sigma_i\lambda_{X_i}X_i``) that becomes ammonia or nitrate depends on the equilibrium of chemical equations dependent on the bottom water ``NO_3``, ``NH_4``, and ``O_2`` concentrations, as well as the total remineralisation and mean degradation rate. The mean first order degredation rate is given by:

```math
k = \frac{\lambda_\text{fast} C_\text{fast} + \lambda_\text{slow} C_\text{slow}}{C_\text{fast} + C_\text{slow}},
```

and, based on [Soetaert2000](@citet), the fraction nitrified (i.e becoming nitrate rather than ammonia) is given by:

```math
\ln\left(C_\text{min}p_{nit}\right) = n_A + n_B\ln C_\text{min}\ln O_2 + n_C * \ln C_\text{min} ^ 2 + n_D * \ln k \ln NH_4 + n_E \ln C_\text{min} + n_F \ln C_\text{min} \ln NH_4.
```

Therefore, the efflux of nitrate and ammonia are given by:

```math
\frac{\partial NO_3}{\partial t} =\frac{N_\text{min}p_\text{nit}}{\Delta z},
```
```math
\frac{\partial NH_4}{\partial t} = \frac{N_\text{min}(1 - p_\text{nit})}{\Delta z},
```

where ``\Delta z`` is the depth of the bottom cell (since ``X_i`` is a surface concentration). This mineralisation also consumes oxygen at a rate:

```math
\frac{\partial O_2}{\partial t} =\frac{C_\text{min}(1 - p_\text{anox}p_\text{solid deposition}) + N_\text{min}p_\text{nit} O:N}{\Delta z},
```

where the constants are given by:

```math
p_\text{solid deposition} = s_A w ^{s_B},
```

with, ``w = s_CD^{s_D}`` where ``D`` is the water depth.

```math
\ln\left(C_\text{min}p_\text{anox}\right) = a_A + a_B\ln C_\text{min} + a_C \ln C_\text{min} ^ 2 + a_D \ln k + a_E \ln O_2 \ln k + a_F \ln NO_3 ^2.
```

The original model of [Soetaert2000](@citet) also includes denitrification terms whereby nitrogen is returned to the water column as dissolved ``N_2``, but we currently do not account for this in order to conserve the nitrogen budget.

### Parameter variable names

| Symbol                  | Variable name              | Units |
|-------------------------|----------------------------|-------|
| ``\lambda_\text{fast}`` | `fast_decay_rate`          | 1 / s |
| ``\lambda_\text{slow}`` | `slow_decay_rate`          | 1 / s |
| ``f_\text{fast}``       | `fast_fraction`            | -     |
| ``f_\text{slow}``       | `slow_fraction`            | -     |
| ``f_\text{ref}``        | `refactory_fraction`       | -     |
| ``n_i``                 | `nitrate_oxidation_params` | -     |
| ``a_i``                 | `anoxic_param`             | -     |
| ``s_i``                 | `solid_dep_params`         | -     |

All parameters are given in [Parameters](@ref parameters).

## Model conservations

Nitrogen and carbon is conserved between the model domain and sediment, any nitrogen or carbon not returned to the bottom cell is stored in a sediment field.

## Model compatibility

This model is currently only compatible with the [LOBSTER](@ref LOBSTER) biogeochemical model.

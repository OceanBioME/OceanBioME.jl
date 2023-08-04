# [Simple multi-G](@id multi-g)

This model, proposed by [Soetaert2000](@citet), is a "G class" model that evolves carbon and nitrogen in three classes (fast, slow and refectory). The model is also only compatible with the LOBSTER biogeochemical model with carbonate chemistry, oxygen, and variable redfield options on. You also must ensure that the `open_bottom` option is on for particles to leave the bottom of the domain to the sediment model.

It is straightforward to set up:

```jldoctest simplemultig; filter = r".*@ OceanBioME.Boundaries.Sediments.*"
using OceanBioME, Oceananigans, OceanBioME.Sediments

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200))

sediment_model = SimpleMultiG(; grid)

# output
┌ Warning: Sediment models are an experimental feature and have not yet been validated.
└ @ OceanBioME.Boundaries.Sediments ~/OceanBioME.jl/src/Boundaries/Sediments/simple_multi_G.jl:102
[ Info: This sediment model is currently only compatible with models providing NH₄, NO₃, O₂, and DIC.
Single-layer multi-G sediment model (Float64)
```

You may optionally specify the model parameters. This can then be passed in the setup of a BGC model:

```julia
biogeochemistry = LOBSTER(; grid, 
                            carbonates = true, oxygen = true, variable_redfield = true, 
                            open_bottom = true, 
                            sediment_model)
```

### Model equations

This model evolved the carbon and nitrogen components of three liability classes: fast, slow, and refractory. Each component is remineralised with first order decay so evolved like:

$\frac{dX_i}{dt} = F_{X_i} - \lambda_iX_i.$

For the fast and slow classes $\lambda$ is a positive, non-zero, rate constant, and for the refractory class it is $0$. $F_{X_i}$ is the flux, and the flux for each class is simply a constant fraction of the total flux:

$F_{X_i} = f_iF_X.$

The fraction of remineralised sediment ($X_\text{min} = \Sigma_i\lambda_{X_i}X_i$) that becomes ammonia or nitrate depends on the equilibrium of chemical equations dependent on the bottom water $NO_3$, $NH_4$, and $O_2$ concentrations, as well as the total remineralisation and mean degradation rate. The mean first order degredation rate is given by:

$k = \frac{\lambda_\text{fast} C_\text{fast} + \lambda_\text{slow} C_\text{slow}}{C_\text{fast} + C_\text{slow}},$

and, based on [Soetaert2000](@citet), the fraction nitrified (i.e becoming nitrate rather than ammonia) is given by:

$\ln\left(C_\text{min}p_{nit}\right) = A + B\ln C_\text{min}\ln O_2 + C * \ln C_\text{min} ^ 2 + D * \ln k \ln NH_4 + E \ln C_\text{min} + F \ln C_\text{min} \ln NH_4.$

Therefore, the efflux of nitrate and ammonia are given by:

$\frac{\partial NO_3}{\partial t} =\frac{N_\text{min}p_\text{nit}}{\Delta z},$
$\frac{\partial NH_4}{\partial t} = \frac{N_\text{min}(1 - p_\text{nit})}{\Delta z},$

where $\Delta z$ is the depth of the bottom cell (since $X_i$ is a surface concentration). This mineralisation also consumes oxygen at a rate:

$\frac{\partial O_2}{\partial t} =\frac{C_\text{min}(1 - p_\text{anox}p_\text{solid deposition}) + N_\text{min}p_\text{nit} O:N}{\Delta z},$

where the constants are given by:

$p_\text{solid deposition} = 0.233 w ^{0.336},$

with, $w = 982D^{-1.548}$ where $D$ is the water depth.

$\ln\left(C_\text{min}p_\text{anox}) = A + B\ln C_\text{min} + C \ln C_\text{min} ^ 2 + D \ln k + E \ln O_2 \ln k + F \ln NO_3 ^2.$

All parameters are given in [Parameters](@ref parameters). The original model of [Soetaert2000](@citet) also includes denitrification terms whereby nitrogen is returned to the water column as dissolved $N_2$, but we currently do not account for this in order to conserve the nitrogen budget.

### Model conservations

Nitrogen and carbon is conserved between the model domain and sediment, any nitrogen or carbon not returned to the bottom cell is stored in a sediment field.

### Model compatibility

This model is currently only compatible with the [LOBSTER](@id LOBSTER) biogeochemical model.
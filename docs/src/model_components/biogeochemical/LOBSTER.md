# [The Lodyc-DAMTP Ocean Biogeochemical Simulation Tools for Ecosystem and Resources (LOBSTER) model](@id LOBSTER)

LOBSTER is a medium complexity BGC model with seven core prognostic variables: phytoplankton, zooplankton, small and large detritus, nitrates, ammonia, and dissolved organic matter. LOBSTER was originally proposed by [Levy2005](@citet) and subsequently added to by [Levy2001](@citet), [Resplandy2009](@citet), [Karleskind2011](@citet), and [Resplandy2012](@citet).

Additionally, this implementation of LOBSTER optionally models simple carbonate chemistry (`DIC` and `Alk`alinity), `Oxy`gen, and variable redfield ratios for the now dissolved and particulate organic groups (which then allows carbon to be conserved). For details see [StrongWrightInPrep](@citet). These are activated in the model setup, for example:

```jldoctest
julia> using OceanBioME, Oceananigans

julia> grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

julia> bgc_model = LOBSTER(; grid, carbonates = true)
Lodyc-DAMTP Ocean Biogeochemical Simulation Tools for Ecosystem and Resources (LOBSTER) model (Float64) 
 Light Attenuation Model: 
    └── Two-band light attenuation model (Float64)
 Optional components:
    ├── Carbonates ✅ 
    ├── Oxygen ❌ 
    └── Variable Redfield Ratio ❌
 Sinking Velocities:
    ├── sPOM: 0.0 to -3.47e-5 m/s 
    └── bPOM: 0.0 to -0.0023148148148148147 m/s
```

# [Air-sea gas exchange](@id air-sea-gas)

We currently have one air-sea gas exchange model implemented. The model, proposed by [Wanninkhof1992](@citet), calculates the solubility of the gas in the water dependent on the temperature and salinity, and calculates the flux depending on the solubility and mixing from the wind.

Currently, the parameters for CO₂ and oxygen are included, but it would be very straightforward to add the parameters given in the original publication for other gases (e.g. inert tracers of other nutrients such as N₂). We also currently have a very simple formulation of the gas transfer velocity which depends on an average wind speed, but it would straightforwardly be expanded to permit variable wind speed (e.g. to simulate enhanced exchange from storms).

It is straightforward to set up a boundary as an air-sea gas exchange:

```jldoctest gasexchange
julia> using OceanBioME

julia> CO₂_flux = GasExchange(; gas = :CO₂)
FluxBoundaryCondition: ContinuousBoundaryFunction gasexchange_function at (Nothing, Nothing, Nothing)
```

Where the symbol specifies the exchanged gas (currently `:CO₂` or `:O₂`). This can then be passed in the setup of a BGC model, for example:

```jldoctest gasexchange
julia> using Oceananigans

julia> model = NonhydrostaticModel(; grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200)),
                                     biogeochemistry = LOBSTER(; grid, carbonates = true),
                                     boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ))
NonhydrostaticModel{CPU, RectilinearGrid}(time = 0 seconds, iteration = 0)
├── grid: 3×3×30 RectilinearGrid{Float64, Periodic, Periodic, Bounded} on CPU with 3×3×3 halo
├── timestepper: QuasiAdamsBashforth2TimeStepper
├── tracers: (NO₃, NH₄, P, Z, sPOM, bPOM, DOM, DIC, Alk)
├── closure: Nothing
├── buoyancy: Nothing
└── coriolis: Nothing
```

If the temperature and salinity are not included in the model they can be passed as functions
(or even anonymous functions):

```jldoctest gasexchange
julia> T_function(x, y, z, t) = 12.0
T_function (generic function with 1 method)

julia> CO₂_flux = GasExchange(; gas = :CO₂, temperature = T_function, salinity = (args...) -> 35)
FluxBoundaryCondition: ContinuousBoundaryFunction gasexchange_function at (Nothing, Nothing, Nothing)
```

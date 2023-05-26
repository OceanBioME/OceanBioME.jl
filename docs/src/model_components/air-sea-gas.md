# [Air-sea gas exchange](@id air-sea-gas)

We currently have one air-sea gas exchange model implemented. The model, proposed by [Wanninkhof1992](@cite), calculates the solubility of the gas in the water dependent on the temperature and salinity, and calculates the flux depending on the solubility and mixing from the wind.

Currently, the parameters for CO₂ and oxygen are included, but it would be very straightforward to add the parameters given in the original publication for other gases (e.g. inert tracers of other nutrients such as N₂). We also currently have a very simple formulation of the gas transfer velocity which depends on an average wind speed, but it would straightforwardly be expanded to permit variable wind speed (e.g. to simulate enhanced exchange from storms).

It is straightforward to set up a boundary as an air-sea gas exchange:

```julia
CO₂_flux = GasExchange(; gas = :CO₂)
```

Where the symbol specifies the exchanged gas (currently `:CO₂` or `:O₂`). This can then be passed in the setup of a BGC model, for example:

```julia
model = NonhydrostaticModel(; grid,
                              biogeochemistry = LOBSTER(; grid,
                                                          carbonates = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ),)
```

If the temperature and salinity are not included in the model they can be passed as functions:
```julia
CO₂_flux = GasExchange(; gas = :CO₂, temperature = t_function, salinity = s_function)
```
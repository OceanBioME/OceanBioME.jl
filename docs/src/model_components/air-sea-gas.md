# [Air-sea gas exchange](@id air-sea-gas)

We currently have one air-sea gas exchange model implemented. The model, proposed by [Wanninkhof1992](@citet), calculates the solubility of the gas in the water dependent on the temperature and salinity, and calculates the flux depending on the solubility and mixing from the wind.

Currently, the parameters for CO₂ and oxygen are included, but it would be very straightforward to add the parameters given in the original publication for other gases (e.g. inert tracers of other nutrients such as N₂). We also currently have a very simple formulation of the gas transfer velocity which depends on an average wind speed, but it would straightforwardly be expanded to permit variable wind speed (e.g. to simulate enhanced exchange from storms).

It is straightforward to set up a boundary as an air-sea gas exchange:

```@setup gasexchange
using OceanBioME
CO₂_flux = GasExchange(; gas = :CO₂)
using Oceananigans

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

model = NonhydrostaticModel(; grid,
                              biogeochemistry = LOBSTER(; grid, carbonates = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ),
                              tracers = (:T, :S))
```


```@example gasexchange
using OceanBioME

CO₂_flux = GasExchange(; gas = :CO₂)
```

Where the symbol specifies the exchanged gas (currently `:CO₂` or `:O₂`). This can then be passed in the setup of a BGC model, for example:

```@example gasexchange
using Oceananigans

grid = RectilinearGrid(size=(3, 3, 30), extent=(10, 10, 200));

model = NonhydrostaticModel(; grid,
                              biogeochemistry = LOBSTER(; grid, carbonates = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ),
                              tracers = (:T, :S))
```

## Model equations

The gas flux is given by:

```math
F = k(C_w - \alpha C_a),
```

where ``C_w`` is the concentration in the water, ``C_a`` the concentration in the air, ``\alpha`` the Oswald solubility coefficient, and ``k`` the gas transfer velocity. For carbon dioxide the flux is modified to:

```math
F = k\beta\rho_o(pCO_{2}(water) - pCO_{2}(air)),
```

where ``pCO_{2}(water)`` and ``pCO_{2}(air)`` are the partial pressure of carbon dioxide in the water and air, ``\beta`` is the Bunsen Solubility Coefficient, and ``\rho_o`` is the density of the water.
``pCO_{2}(water)`` is diagnosed from the dissolved inorganic carbon and alkalinity of the water using the [carbon chemistry](@ref carbon-chemistry) model.

The gas transfer velocity is parameterised by the wind speed and Schmidt number, which in turn is parameterised by the temperature and salinity. The gas transfer velocity is given by:

```math
k = 1.08\times10^{-6}u^2\left(\frac{Sc}{660}\right)^{-1/2},
```

where ``u`` is the winds speed 10m above the surface, and Sc is the Schmidt number parameterised as:

```math
Sc = A - BT + CT^2 - DT^3,
```

where ``T`` is temperature in Kelvin and the other parameters are dependent on the gas type and given in [Parameters](@ref parameters).

The solubilities are given by:

```math
\alpha = 0.00367 T \exp{A_1 + 100\frac{A_2}{T} + A_3 \ln{\frac{T}{100}} + S\left(B_1 + \frac{B_2}{T} + \frac{B_3}{T^2}\right)},
```

and

```math
\beta = \exp{A_1 + 100\frac{A_2}{T} + A_3 \ln{\frac{T}{100}} + S\left(B_1 + \frac{B_2}{T} + \frac{B_3}{T^2}\right)},
```

where ``S`` is salinity in practical units and the other default parameters are given in [Parameters](@ref parameters).

# [Biogeochemical Models](@id bgc_models)

Biogeochemical (BGC) models are generally implemented as [forced tracers in Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/model_setup/tracers/), or can be run as box models. All BGC models should be setup in the same way so that they can easily be substituted for each other.

For details of the BGC models currently implemented please see the following pages.

## Oceananigans setup
At the simplest level all that is required to setup a BGC model is the grid the model will run on, which we pass to a setup function:
```
bgc = Setup.Oceananigans(MODEL_NAME, grid, MODEL_PARAMETERS) 
```
Here you must specify the `MODEL_NAME`, and parameters. The default parameters can be accessed at `MODEL_NAME.defaults`.

You may also wish to specify additional boundary conditions for any of the BGC model tracers (which you can see at `MODEL_NAME.tracers`). A simple example is some addition of nutrients from the bottom:
```
NO₃_upwelling = FluxBoundaryCondition(0.1)
bgc = Setup.Oceananigans(MODEL_NAME, grid, MODEL_PARAMETERS, bottom_boundaries=(NO₃ = NO₃_upwelling,)) 
```

More likely you will want to use [air-sea exchange](@ref air-sea-gas) at the top or [sediment](@ref sediments) at the bottom.

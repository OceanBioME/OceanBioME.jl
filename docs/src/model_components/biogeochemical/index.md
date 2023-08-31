# [Biogeochemical Models](@id bgc_models)

Biogeochemical (BGC) models can be used within the [Oceananigans biogeochemistry framework](https://github.com/CliMA/Oceananigans.jl/pull/2802) or as stand alone box models. All BGC models should be setup in the same way so that they can easily be substituted for each other. You can easily implement a different model (or a variation on a current model) by following the guide [here](@ref model_implementation).

For details of the BGC models currently implemented please see the following pages.

## Oceananigans setup
At the simplest level all that is required to setup an existing OceanBioME BGC model is to pass it to the Oceananigans model setup:
```julia
model = NonhydrostaticModel(; grid,
                              ...,
                              biogeochemistry = MODEL_NAME(; grid))
```
Where `MODEL_NAME` is the name of the model, you may also need to pass other parameters like:
```julia
MODEL_NAME(; grid, growth_rate = 10.0)
```

This will set up the required tracers and auxiliary fields, and you may also set boundary conditions or additional forcing through the normal Oceananigans setup. 

Models usually have a default [light attenuation model](@ref light) specified, these may be substituted easily by passing different models as parameters as above.

Our models are implemented in an abstract framework `Biogeochemistry` which contains `underlying_biogeochemistry`, `light_attenuation`, `sediment`, and `modifiers`.
This is automatically setup for existing BGC models, but may also be used to couple any BGC model with light attenuation and sediments.
See the [implementation page](@ref model_implementation) for some more information on how to couple other models.
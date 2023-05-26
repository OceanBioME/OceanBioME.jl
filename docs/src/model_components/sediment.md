# [Sediment](@id sediment)

We currently have one single layer sediment model implemented. The model, proposed by [Soetaert2000](@cite), is a "G class" model that evolves carbon and nitrogen in three classes (fast, slow and refectory). The model is also only compatible with the LOBSTER biogeochemical model with carbonate chemistry, oxygen, and variable redfield options on. You also must ensure that the `open_bottom` option is on for particles to leave the bottom of the domain to the sediment model.

It is straightforward to set up:

```julia
sediment_model = SimpleMultiG(grid)
```

You may optionally specify the model parameters. This can then be passed in the setup of a BGC model:

```julia
biogeochemistry = LOBSTER(; grid, 
                            carbonates = true, oxygen = true, variable_redfield = true, 
                            open_bottom = true, 
                            sediment_model)
```
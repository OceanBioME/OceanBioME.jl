# [Soetaert et al. 2000](@id soetaert)

This is a simple - single layer - sediment model with three classes of particulate matter which decay at different rates (by default fast, slow, and refractory), first proposed in [Soetaert2000](@cite). The default way to set it up is simply:
```
sediment=Boundaries.Sediments.Soetaert.setupsediment(grid)
```
But you may pass optional parameters:
- `POM_fields` - tuple of sinking particulate matter tracer fields,
- `POW_w` - sinking speed of the tracer fields

And others.

As with [air sea flux](@ref air-sea-gas) you then need to tell the BGC model about this, e.g.:
```
bgc = Setup.Oceananigans(:LOBSTER, grid, bottomboundaries=sediment.boundary_conditions)
```
And then add the correct forcing and auxiliary fields to the model:
```
model = NonhydrostaticModel(
...
    forcing =  merge(bgc.forcing, sediment.forcing),
    boundary_conditions = bgc.boundary_conditions,
    auxiliary_fields = merge((PAR=PAR, ), sediment.auxiliary_fields)
)
```
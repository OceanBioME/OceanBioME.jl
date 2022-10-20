# [Light attenuation models](@id light)

Nearly all BGC models require some model of the attenuation of PAR through the water. Usually this depends on the concentration of chlorophyll in the water (in phytoplankton), and may depend on the concentration of colored dissolved organic matter or particulates.

We currently have two models of light attenuation, a two band model by [Karleskind2011](@cite) and the more widely used three band model by [Morel1988](@cite). As the light level is diagnostic of the phytoplankton concentration these models are implemented with the light level as various auxiliary fields which are updated with callbacks. Different BGC models expect different PAR fields (which you will be warned about when the model is setup, or you can find from `MODEL_NAME.required_fields`), you need to specify these, e.g.:
```
PAR_field = Oceananigans.Fields.Field{Center, Center, Center}(grid)
```
and pass them as auxiliary fields to the model:
```
model = NonhydrostaticModel(
...
    auxiliary_fields = (PAR=PAR_field, )
)
```
Finally, once the simulation is defined you need to add the callback to update the field:
```
simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(params, Light.twoBands.defaults, (;surface_PAR)));
```
Where `surface_PAR` is some function `surface_PAR(t)` which gives the light level at the surface at time `t`. 

For the Morel model the setup would instead be:
```
PAR¹, PAR², PAR³ = Field{Center, Center, Center}(grid), Field{Center, Center, Center}(grid), Field{Center, Center, Center}(grid)
...
model = NonhydrostaticModel(
...
    auxiliary_fields = (; PAR¹, PAR², PAR³)
)
...
simulation.callbacks[:update_PAR] = Callback(Light.Morel.update, IterationInterval(1), merge(Light.Morel.defaults, (; PAR⁰)))
```

Since this clearly wouldn't have the correct fields for a model expecting a single PAR field (e.g. LOBSTER) you could also define a custom field and callback which converts them, making the setup:
```
PAR¹, PAR², PAR³, PAR = Field{Center, Center, Center}(grid), Field{Center, Center, Center}(grid), Field{Center, Center, Center}(grid), Field{Center, Center, Center}(grid)
...
model = NonhydrostaticModel(
...
    auxiliary_fields = (; PAR¹, PAR², PAR³, PAR)
)
...
function update_PAR!(sim, params)
    Light.Morel.update(sim, params)
    PAR .= PAR¹ .+ PAR² .+ PAR³
end
simulation.callbacks[:update_PAR] = Callback(update_PAR!, IterationInterval(1), merge(Light.Morel.defaults, (; PAR⁰)))
```
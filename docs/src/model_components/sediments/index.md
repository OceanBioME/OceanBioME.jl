# [Sediment](@id sediment)

Sediment models can be added to biogeochemical models. For details of the BGC models currently implemented please see the following pages. 

Sediment models are usually added by setting up the model and then passing it to the biogeochemical model, for example:
```
sediment = SEDIMENT_MODEL_NAME(; grid)

biogeochemistry = BIOGEOCHEMICAL_MODEL_NAME(; name, sediment, ...)
```

where `SEDIMENT_MODEL_NAME` is the chosen sediment model, and `BIOGEOCHEMICAL_MODEL_NAME` is the chosen biogeochemical model and `...` replaces the other parameters you may wish to pass to the model.

Please note that not all sediment models are compatible with all biogeochemical models. This will be noted on the documentation page for each model.
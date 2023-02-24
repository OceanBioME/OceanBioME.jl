# [Light attenuation models](@id light)

Nearly all BGC models require some model of the attenuation of PAR through the water. Usually this depends on the concentration of chlorophyll in the water (in phytoplankton), and may depend on the concentration of coloured dissolved organic matter or particulates.

We currently have two models of light attenuation, a two band model by [Karleskind2011](@cite) and the more widely used three band model by [Morel1988](@cite). As the light level is diagnostic of the phytoplankton concentration these models are implemented with the light level as various auxiliary fields which are updated with callbacks within the biogeochemical model.

Models requiring light attenuation models will set these up automatically, for example [LOBSTER](@ref LOBSTER) sets `light_attenuation_model = TwoBandPhotosyntheticallyActiveRatiation()`. You may choose others. Additionally, you can pass the surface PAR as a function of horizontal position and time. The default for LOBSTER is `(x, y, t) -> 100*max(0.0, cos(t*π/(12hours)))`.
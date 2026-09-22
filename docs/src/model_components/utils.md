# [Utilities](@id utils)

We provide some utilities that may be useful.

## Time step adaptation
We have added a few additional utilities which extend the capabilities of Oceananigans' time step wizard. For column models where there is no water velocity we have added functions to calculate the advection timescale from the biogeochemical model defined sinking velocities. This could be used by:
```julia
wizard = TimeStepWizard(cfl = 0.2, diffusive_cfl = 0.2, max_change = 2.0, min_change = 0.5, cell_advection_timescale = column_advection_timescale)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
```
Finally, sinking may be more limiting than the normal advective CFL conditions so, we have an additional cell advection timescale defined for 3D models:
```julia
wizard = TimeStepWizard(cfl = 0.6, diffusive_cfl = 0.5, max_change = 1.5, min_change = 0., cell_advection_timescale = sinking_advection_timescale)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
```

## Negative tracer treatment
Negative tracer values can be handled through the `negative_tracers` keyword of `Biogeochemistry` and model constructors. Three treatments are available.

`ScaleNegativeTracers` rescales a conserved tracer group so none of its members remain negative:
```julia
negativity_protection = ScaleNegativeTracers((:P, :Z, :N))
biogeochemistry = Biogeochemistry(...; negative_tracers = negativity_protection)
```
A scale factor can be provided for each component, for example when tracers use different Redfield ratios:
```julia
negativity_protection = ScaleNegativeTracers((:P, :Z, :N); scalefactors = (1, 1, 2))
biogeochemistry = Biogeochemistry(...; negative_tracers = negativity_protection)
```
The conserved group should be chosen to match the model's elemental bookkeeping. For supported model constructors, `ScaleNegativeTracers()` can infer the conserved group or groups from the underlying biogeochemistry.

`ClipNegativeTracers` directly clips negative prognostic tracer values to zero:
```julia
biogeochemistry = Biogeochemistry(...; negative_tracers = ClipNegativeTracers())
```
It optionally accepts `exclude`, a tuple of tracer names (`Symbol`s) that may remain negative.

`IgnoreNegativeTracerValues` leaves the prognostic tracer fields unchanged while biogeochemical processes evaluate negative concentration-tracer values as zero:
```julia
biogeochemistry = Biogeochemistry(...; negative_tracers = IgnoreNegativeTracerValues())
```
Signed environmental tracers such as temperature and salinity remain unchanged. The same treatment is applied to chlorophyll values used by light attenuation.

Treatments can be composed in a tuple when evaluation-time and state-update behaviour are both desired:
```julia
biogeochemistry = Biogeochemistry(...;
    negative_tracers = (IgnoreNegativeTracerValues(), ScaleNegativeTracers((:P, :Z, :N))))
```

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

## Negative tracer detection
As a temporary measure we have implemented a callback to either detect negative tracers and either scale a conserved group, force them back to zero, or throw an error. Please see the numerical implementations' page for details. This can be set up by:
```julia
negativity_protection = ScaleNegativeTracers((:P, :Z, :N))
biogeochemistry = Biogeochemistry(...; modifiers = negativity_protection)
```
You may also pass a scale factor for each component (e.g. in case they have different redfield ratios):
```julia
negativity_protection = ScaleNegativeTracers((:P, :Z, :N); scalefactors = (1, 1, 2))
biogeochemistry = Biogeochemistry(...; modifiers = negativity_protection)
```
Here you should carefully consider which tracers form a conserved group (if at all). If some tracers are in more than one group (e.g. plankton contain both nitrogen and carbon), give all of the groups to the same `ScaleNegativeTracers` so that they are all conserved:
```julia
negativity_protection = ScaleNegativeTracers((nitrogen = (N = 1, P = 1, Z = 1),
                                              carbon = (DIC = 1, P = 6.56, Z = 6.56)))
biogeochemistry = Biogeochemistry(...; modifiers = negativity_protection)
```
Separate `ScaleNegativeTracers` modifiers are applied one after the other, so each would change the totals of the others. The `scale_negatives = true` option of the models does this for the groups given by `OceanBioME.conserved_tracers(biogeochemistry)`. Alternatively, force to zero by:
```julia
negativity_protection = ZeroNegativeTracers()
biogeochemistry = Biogeochemistry(...; modifiers = negativity_protection)
```
The latter optionally takes a named tuple of parameters that may include `exclude`, which can be a tuple of tracer names (Symbols) which are allowed to be negative.

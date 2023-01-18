# [Utilities](@id utils)

We provide some utilities that may be useful.

## Time step adaptation
To automatically adapt the time step length you may add a callback like:
```
simulation.callbacks[:timestep] = Callback(update_timestep!, IterationInterval(1), (w=200/day, c_diff = 0.45, c_adv = 0.45, relaxation=0.75))
```
`c_diff` and `c_adv` are the diffusive and advective [Courant Numbers](https://www.wikiwand.com/en/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition), and the relaxation damps the change in step length as the new step is ``\Delta t_{i+1} = \Delta t_i \left(\frac{C_{max}}{C}\right)^{relaxation}``.

Optionally you can also specify a maximum time step length `Δt_max` and experimentally `c_forcing`, although we do not define this in the usual Courant number way but instead as `max(G/C)`. This reduces the time step when the tendency becomes large/the concentration becomes small in an attempt to prevent numerical error taking tracers below zero, but there is no mathematical reason for the instability to scale like this. We [plan](https://github.com/orgs/OceanBioME/projects/4) on implementing a positivity preserving time stepper in the future which would overcome this issue.

## Negative tracer detection
As a temporary measure we have implemented a callback to either detect negative tracers and either scale a conserved group, force them back to zero, or throw an error. Please see the numerical implementations page for details. This can be set up by:
```
negativity_protection! = ScaleNegativeTracers(tracers = (:P, :Z, :N))
simulation.callbacks[:neg] = Callback(negativity_protection!; callsite = UpdateStateCallsite())
```
You may also pass a scale factor for each component (e.g. in case they have different redfield raitos):
```
negativity_protection! = ScaleNegativeTracers(tracers = (:P, :Z, :D), scalefactors = (P = 1, Z = 1, D = 2))
simulation.callbacks[:neg] = Callback(negativity_protection!; callsite = UpdateStateCallsite())
```
Here you should carefully consider which tracers form a conserved group (if at all). Alternatively, force to zero by:
```
simulation.callbacks[:neg] = Callback(OceanBioME.no_negative_tracers!, callsite = UpdateStateCallsite())
```
or throw an error:
```
simulation.callbacks[:neg] = Callback(OceanBioME.error_on_neg!, callsite = UpdateStateCallsite())
```
The latter two both optionally take a named tuple of parameters which may include `exclude` which can be a tuple of tracer names (Symbols) which are allowed to be negative.
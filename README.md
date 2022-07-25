# BGC.jl (placeholder name, changed to reduce confusion)

## To use (before release):

- Activate julia with the `--project` flag pointed to this directory, i.e. `julia --project`
- Ensure all the dependencies are installed by typing `] instantiate` or `using Pkg; Pkg,instantiate()` (if the dependencies are not currently listed properly you may have to manually add them)
- Now you can use the package with `using BGC` as usual

## Usage overview:

1. Define the physical dependencies: `T`, `S`, and `PAR`. Currently, these can be functions of `(x, y, z, t)`, or, for `T` and `S` tracer fields (performs much slower). Additionally, `PAR` can be an auxiliary field calculated by callbacks provided by BGC.jl. For example:
```
t_function(x, y, z, t) = temperature_itp(mod(t, 364days)) .+ 273.15
s_function(x, y, z, t) = salinity_itp(mod(t, 364days))
surface_PAR(t) = surface_PAR_itp(mod(t, 364days))

PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)
```
This latter example requires that the grid is already defined.

2. Pass these fields and the grid to a setup function (currently Lobster is the only one implimented)
```
bgc = Lobster.setup(grid, params, (T=t_function, S=s_function, PAR=PAR))
```
This returns a model with fields `tracers`, `forcing`, and `boundary_conditions`.

3. Pass the result to an Oceananigans model (along with the auxiliary fields if necessary)
```
model = NonhydrostaticModel(...
                            tracers = (:b, bgc.tracers...),
                            ...
                            forcing =  bgc.forcing,
                            boundary_conditions = merge((u=u_bcs, b=buoyancy_bcs), bgc.boundary_conditions),
                            auxiliary_fields = (PAR=PAR, ))
```

4. Finally, once you have otherwise defined the simulation, setup the callback to update the PAR (here Light.update_2λ! is a function provided by BGC.jl that requires a surface PAR interpolant be passed to it. (At some point this will need to be updated to be an x, y, t interpolation).
```
simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)))

```

## Notes on library structure
I think having the models (e.g. LOBSTER, NPZ, etc.) as submodules which are exported by the library works quite well because it allows their forcing functions etc. to be accessed/extended with no ambiguity about which model they belong to. For example, if I want to know about the phytoplankton forcing in the Lobster model I can type:
```
?Lobster.P_forcing
  No documentation found.

  BGC.Lobster.P_forcing is a Function.

  # 2 methods for generic function "P_forcing":
  [1] P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) in BGC.Lobster at /home/jago/Documents/Projects/Lobster/src/Lobster.jl:21
  [2] P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) in BGC.Lobster at /home/jago/Documents/Projects/Lobster/src/Lobster.jl:34

```
And if other models existed this could be e.g. `NPZ.P_forcing`. 

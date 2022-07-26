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

## Particles
To include particles which interact with the BGC model setup the model as above, but before the Oceananigans model is defined set up the particles by:

1. Define a function of `(x, y, z, t, property dependencies..., params, Δt)`, for example:
```
function sugarkelpequations(x, y, z, t, A, N, C, NO₃, irr, params, Δt)
  dA, dN, dC, j = SugarKelp.equations(A, N, C, params.urel, params.temp(mod(t, 364days), 
                                                irr, NO₃, params.λ[1+floor(Int, mod(t, 364days)/day)], params.resp_model, Δt, params.paramset)
  return (A = dA / (60*60*24), N = dN / (60*60*24), C = dC / (60*60*24), j = (A+Δt*dA / (60*60*24)) * j / (60*60*24*14*0.001))#fix units
end
lat=57.5
λ_arr=SugarKelp.gen_λ(lat)
```
2. Next create the particle structure as is done with Oceananigans
```
struct Kelp
  #position
  x :: AbstractFloat
  y :: AbstractFloat
  z :: AbstractFloat
  #velocity
  u :: AbstractFloat
  v :: AbstractFloat
  w :: AbstractFloat
  #properties
  A :: AbstractFloat
  N :: AbstractFloat
  C :: AbstractFloat
  j :: AbstractFloat
  #tracked fields
  NO₃  :: AbstractFloat
  PAR :: AbstractFloat
end
... #define initial value arrays
particles = StructArray{Kelp}((x₀ₖ, y₀ₖ, z₀ₖ, u₀ₖ, v₀ₖ, w₀ₖ, a₀, n₀, c₀, j₀, NO₃₀, PAR₀))
```
3. Define the "source fields" which are those tracked by the particle (usually arguments of the function), and "sink fields" which are the connection between particle properties and tracers
```
source_fields = ((tracer=:NO₃, property=:NO₃, scalefactor=1.0), (tracer=:PAR, property=:PAR, scalefactor=1.0))
sink_fields = ((tracer=:NO₃, property=:j, scalefactor=-1.0), )
``` 
4. Call the particle setup function to define Oceananigans LagrangianParticles
``` 
kelp_particles = Particles.setup(particles, sugarkelpequations, 
                                        (:A, :N, :C, :NO₃, :PAR), #forcing function property dependencies
                                        (urel=0.15, temp=temperature_itp, λ=λ_arr, resp_model=2, paramset=SugarKelp.broch2013params), #forcing function parameters
                                        (:A, :N, :C), #forcing function integrals
                                        (:j, ), #forcing function diagnostic fields
                                        source_fields,
                                        sink_fields,
                                        100.0)#density
``` 
5. Finally, setup the model with the particles
``` 
model = NonhydrostaticModel(advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            grid = grd,
                            tracers = (:b, bgc.tracers...),
                            coriolis = FPlane(f=1e-4),
                            buoyancy = BuoyancyTracer(), 
                            closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                            forcing =  bgc.forcing,
                            boundary_conditions = merge((u=u_bcs, b=buoyancy_bcs), bgc.boundary_conditions),
                            auxiliary_fields = (PAR=PAR, ),
                            particles = kelp_particles)
``` 
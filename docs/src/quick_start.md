# Quick start
This code will run one month of a single column, 7 variable (P, Z, D, DD, DOM, NO₃, NH₄) biogeochemical situation with constant forcing.

```@meta
DocTestSetup = quote
    using OceanBioME, Oceananigans, Plots, NetCDF
    using Oceananigans.Units: days, minute
end
```

``` jldoctest quickstart
grid = RectilinearGrid(size=(1, 1, 10), extent=(1, 1, 200), topology=(Periodic, Periodic, Bounded))

bgc = Setup.Oceananigans(:LOBSTER, grid, LOBSTER.defaults) 

PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

model = NonhydrostaticModel(
    grid = grid,
    tracers = bgc.tracers,
    forcing = bgc.forcing,
    boundary_conditions = bgc.boundary_conditions,
    auxiliary_fields = (; PAR)
)

set!(model, P = 0.001, Z = 0.001, NO₃ = 1, NH₄ = 0.01)
surface_PAR(t) = 10

simulation = Simulation(model; Δt=1minute, stop_time=30days)
simulation.callbacks[:update_PAR] = Callback(Light.twoBands.update!, IterationInterval(1), merge(LOBSTER.defaults, Light.twoBands.defaults, (;surface_PAR)))
simulation.output_writers[:profiles] = NetCDFOutputWriter(model, model.tracers[bgc.tracers], filename="quickstart.nc", schedule=TimeInterval(0.5days), overwrite_existing=true)

run!(simulation)

# output
┌ Warning: This model requires (:PAR,) to be separatly defined (as tracer or auxiliary fields)
└ @ OceanBioME.Setup ~/Documents/Projects/Lobster/src/Utils/Setup.jl:124
┌ Warning: This model requires () to be separatly defined in addition to the default parameters (MODEL_NAME.defaults)
└ @ OceanBioME.Setup ~/Documents/Projects/Lobster/src/Utils/Setup.jl:130
[ Info: Initializing simulation...
[ Info:     ... simulation initialization complete (666.735 ms)
[ Info: Executing initial time step...
[ Info:     ... initial time step complete (4.200 ms).
[ Info: Simulation is stopping. Model time 30 days has hit or exceeded simulation stop time 30 days.

```
This isn't quite as simple as it could be as it records the output so that we can visualize it:

``` jldoctest quickstart

times = ncread("quickstart.nc", "time")

phytoplankton = ncread("quickstart.nc", "P")
nitrates = ncread("quickstart.nc", "NO₃")

hm1 = heatmap(times/days, grid.zᵃᵃᶜ[1:10], phytoplankton[1, 1, :, :], xlabel="Day", ylabel="Depth (m)", title="Phytoplankton (mmol N/m³)")
hm2 = heatmap(times/days, grid.zᵃᵃᶜ[1:10], nitrates[1, 1, :, :], xlabel="Day", ylabel="Depth (m)", title="Nitrate (mmol N/m³)")

plot(hm1, hm2)
savefig("docs/src/img/quickstart.png");
# output
```
![Heatmaps of phytoplankton and nitrate concentration showing phytoplankton growth and nitrate consumption](img/quickstart.png)

OceanBioME provides the tools to add to this, for example adding a carbonate chemistry model, or sediment at the bottom of the model. Please have a look at the rest of the examples to explore these options.
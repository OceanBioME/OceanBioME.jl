# Quick start

OceanBioME provides biogeochemical models to plug into [Oceananigans](https://github.com/CliMA/Oceananigans.jl), for example this code will run one month of a single column, 7 variable (P, Z, sPOM, bPOM, DOM, NO₃, NH₄) biogeochemical situation with constant forcing.

First we need to check we have the required dependencies:
```julia
using Pkg
Pkg.add(["OceanBioME", "Oceananigans"])
```

```@example quickstart
using OceanBioME, Oceananigans
using Oceananigans.Units

grid = RectilinearGrid(size = 10, extent = 200meters, topology = (Flat, Flat, Bounded))

model = NonhydrostaticModel(; grid, biogeochemistry = LOBSTER(grid))

set!(model, P = 0.001, Z = 0.001, NO₃ = 1, NH₄ = 0.01)

simulation = Simulation(model, Δt = 1minute, stop_time = 30days)

simulation.output_writers[:profiles] = JLD2Writer(model, model.tracers,
                                                  filename = "quickstart.jld2",
                                                  schedule = TimeInterval(0.5days),
                                                  overwrite_existing = true)

run!(simulation)
```

We can then visualize it, first check the required packages are installed:

```julia
Pkg.add("CairoMakie")
```

and then load the data and plot:

```@example quickstart
using CairoMakie

phytoplankton = FieldTimeSeries("quickstart.jld2", "P")
nitrates = FieldTimeSeries("quickstart.jld2", "NO₃")

_, _, z = nodes(nitrates)

fig = Figure()

axis_kwargs = (xlabel = "Day", ylabel = "Depth (m)")
ax1 = Axis(fig[1, 1]; title = "Phytoplankton (mmol N/m³)", axis_kwargs...)
ax2 = Axis(fig[1, 2]; title = "Nitrate (mmol N/m³)", axis_kwargs...)

hm1 = heatmap!(ax1, phytoplankton.times / day, z, interior(phytoplankton , 1, 1, :, :)')
hm2 = heatmap!(ax2,      nitrates.times / day, z, interior(nitrates, 1, 1, :, :)')

fig
```

OceanBioME provides the tools to add to this, for example adding a carbonate chemistry model, or sediment at the bottom of the model. Please have a look at the rest of the examples to explore these options.

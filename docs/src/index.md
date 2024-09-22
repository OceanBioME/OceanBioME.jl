# *Ocean* *Bio*geochemical *M*odelling *E*nvironment - OceanBioME

OceanBioME.jl is a fast and flexible ocean biogeochemical modelling environment. It is highly modular and is designed to make it easy to implement and use a variety of biogeochemical and physical models. OceanBioME is built to be coupled with physics models from [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) allowing simulations across a wide range of spatial scales ranging from a global hydrostatic free surface model to non-hydrostatic large-eddy simulations. OceanBioME was designed specifically for ocean carbon dioxide removal applications. Notably, it includes active particles which allow individual-based models to be seamlessly coupled with the flow physics, ecosystem models, and carbonate chemistry.

OceanBioME.jl currently provides a core of several biogeochemical models Nutrient--Phytoplankton--Zooplankton--Detritus ([NPZD](@ref NPZD)), [LOBSTER](https://doi.org/10.1029/2004JC002588), a medium complexity model, and an early implementation of [PISCES](https://www.pisces-community.org/), a complex model. It also provides essential utilities like air-sea gas exchange models to provide appropriate top boundary conditions, a carbon chemistry model for computing the pCO₂, and sediment models to for the benthic boundary. 

OceanBioME.jl includes a framework for integrating the growth of biological/active Lagrangian particles which move around and can interact with the (Eulerian) tracer fields - for example, consuming nutrients and carbon dioxide while releasing dissolved organic material. A growth model for sugar kelp is currently implemented using active particles, and this model can be used in a variety of dynamical scenarios including free-floating or bottom-attached particles.

## Quick install

OceanBioME is a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/).

2. Launch Julia and type

```julia
julia> using Pkg
julia> Pkg.add("OceanBioME")
```

!!! compat "Julia 1.10"

    OceanBioME.jl requires Julia version 1.10 or later.

## Running your first model

As a simple example lets run a Nutrient-Phytoplankton-Zooplankton-Detritus (NPZD) model in a two-dimensional simulation of a buoyancy front. This example requires Oceananigans, so we install that first via:

```julia
using Pkg; Pkg.add("Oceananigans")
```

and then:

```@example quickstart
using OceanBioME, Oceananigans
using Oceananigans.Units

grid = RectilinearGrid(CPU(), size = (160, 32), extent = (10000meters, 500meters), topology = (Bounded, Flat, Bounded))

biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid) 

model = NonhydrostaticModel(; grid, biogeochemistry,
                              advection = WENO(; grid),
			                  closure = AnisotropicMinimumDissipation(),
			                  buoyancy = SeawaterBuoyancy(constant_salinity = true))

@inline front(x, z, μ, δ) = μ + δ * tanh((x - 7000 + 4 * z) / 500)

Pᵢ(x, z) = ifelse(z > -50, 0.03, 0.01)
Nᵢ(x, z) = front(x, z, 2.5, -2)
Tᵢ(x, z) = front(x, z, 9, 0.05)

set!(model, N = Nᵢ, P = Pᵢ, Z = Pᵢ, T = Tᵢ)

simulation = Simulation(model; Δt = 50, stop_time = 4days)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "buoyancy_front.jld2",
                                                       schedule = TimeInterval(24minute),
                                                       overwrite_existing = true)

run!(simulation)
```

We can then load the saved output and visualize it:

```@example quickstart
T = FieldTimeSeries("buoyancy_front.jld2", "T")
N = FieldTimeSeries("buoyancy_front.jld2", "N")
P = FieldTimeSeries("buoyancy_front.jld2", "P")

xc, yc, zc = nodes(T)

times = T.times

using CairoMakie

n = Observable(1)

T_lims = (8.94, 9.06)
N_lims = (0, 4.5)
P_lims = (0.007, 0.02)

Tₙ = @lift interior(T[$n], :, 1, :)
Nₙ = @lift interior(N[$n], :, 1, :)
Pₙ = @lift interior(P[$n], :, 1, :)

fig = Figure(size = (1000, 520), fontsize = 20)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 770, yticks = [-400, -200, 0])
ax1 = Axis(fig[1, 1]; title = "Temperature (°C)", axis_kwargs...)
ax2 = Axis(fig[2, 1]; title = "Nutrients concentration (mmol N / m³)",axis_kwargs...)
ax3 = Axis(fig[3, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)

hm1 = heatmap!(ax1, xc, zc, Tₙ, colorrange = T_lims, colormap = Reverse(:lajolla), interpolate = true)
hm2 = heatmap!(ax2, xc, zc, Nₙ, colorrange = N_lims, colormap = Reverse(:bamako), interpolate = true)
hm3 = heatmap!(ax3, xc, zc, Pₙ, colorrange = P_lims, colormap = Reverse(:bamako), interpolate = true)

Colorbar(fig[1, 2], hm1, ticks = [8.95, 9.0, 9.05])
Colorbar(fig[2, 2], hm2, ticks = [0, 2, 4])
Colorbar(fig[3, 2], hm3, ticks = [0.01, 0.02, 0.03])

rowgap!(fig.layout, 0)

record(fig, "buoyancy_front.mp4", 1:length(times)) do i
    n[] = i
end

nothing #hide
```

![buoyancy_front](buoyancy_front.mp4)

In the example above, `OceanBioME.jl` provides the `biogeochemistry` and everything else is taken care of by `Oceananigans.jl`. For comprehensive documentation of the physics modelling see [Oceananigans' Documentation](https://clima.github.io/OceananigansDocumentation/stable/); for biogeochemistry and other features we provide read below.

## Places to find OceanBioME information

* This documentation, which provides
    * documented examples (browse them starting, e.g., from the [single-column model](@ref OneD_column)),
    * explanations of model implementation methods,
    * details of currently implemented models, and
    * a [library](@ref library_api) documenting all user-facing objects and functions.

* [Discussions on the OceanBioME github](https://github.com/OceanBioME/OceanBioME.jl/discussions)
  
    If you've got a question or something to talk about, don't hesitate to [start a new discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions/new?)!

* [Issues](https://github.com/OceanBioME/OceanBioME.jl//issues) and [pull requests](https://github.com/OceanBioME/OceanBioME.jl/pulls) also contain lots of information about problems we've found, solutions we're trying to implement, and ideas for the future.

## Getting in touch

Whether you need help getting started with OceanBioME, found a bug, want OceanBioME to be more expanded, or just want to chat about our project, you can:

* [Start a discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions). 
* [Open an issue](https://github.com/OceanBioME/OceanBioME.jl/issues). Issues are best if you think the OceanBioME source code needs attention: a bug, a sign error, an important missing feature, or a typo in the documentation.

## Citing

If you use OceanBioME as part of your research, teaching, or other activities, we would be grateful if you could cite our work below and mention OceanBioME by name, as well as citing and [acknowledging Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/#Citing) as without them this package would not be possible.

```bibtex
@article{OceanBioMEJOSS,
  doi = {10.21105/joss.05669},
  url = {https://doi.org/10.21105/joss.05669},
  year = {2023},
  publisher = {The Open Journal},
  volume = {8},
  number = {90},
  pages = {5669},
  author = {Jago Strong-Wright and Si Chen and Navid C. Constantinou and Simone Silvestri and Gregory LeClaire Wagner and John R. Taylor},
  title = {{OceanBioME.jl: A flexible environment for modelling the coupled interactions between ocean biogeochemistry and physics}},
  journal = {Journal of Open Source Software}
}
```


## Funding

OceanBioME.jl is supported through grants from the [Center for Climate Repair at Cambridge](https://www.climaterepair.cam.ac.uk/) and the [Gordon and Betty Moore Foundation](https://www.moore.org/).

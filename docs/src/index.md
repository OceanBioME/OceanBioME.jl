# *Ocean* *Bio*geochemical *M*odelling *E*nvironment - OceanBioME

OceanBioME.jl is a fast and flexible ocean biogeochemical modelling environment. It is highly modular and is designed to make it easy to implement and use a varitey of biogeochemical and physical models. OceanBioME is built to be coupled with physics models from [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) allowing simulations across a wide range of spatial scales ranging from a global hydrostatic free surface model to nonhydrostatic large-eddy simulations. OceanBioME was designed specifically for ocean CDR applications. Notably, it includes active particles which allow individual-based models to be seamlessly coupled with the flow physics, ecosystem models, and carbonate chemistry.

OceanBioME.jl is supported through grants from the [Center for Climate Repair at Cambridge](https://www.climaterepair.cam.ac.uk/) and the [Gordon and Betty Moore Foundation](https://www.moore.org/). 

OceanBioME.jl currently provides a core of several biogeochemical models (NPZD and [LOBSTER](https://doi.org/10.1029/2004JC002588), a medium complexity model, and [PISCES](https://doi.org/10.5194/gmd-8-2465-2015) in an early stage of testing), air-sea gas exchange models to provide appropriate top boundary conditions, and sediment models for the benthic boundary (under development).

OceanBioME includes a framework for integrating the growth of biological/active Lagrangian particles which move around and can interact with the (Eulerian) tracer fields - for example, consuming nutrients and carbon dioxide while releasing dissolved organic material. A growth model for sugar kelp is currently implemented using active particles, and this model can be used in a variety of dynamical scenarios including free-floating or bottom-attached particles.

## Quick install

OceanBioME is a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/).

2. Launch Julia and type

```julia
julia> using Pkg
julia> Pkg.add("OceanBioME")
```

## Running your first model
As a simple example let's run a Nutrient-Phytoplankton-Zooplankton-Detritus (NPZD) model in a two-dimensional simulation of a buoyancy front:

```@example quickstart
using OceanBioME, Oceananigans
using Oceananigans.Units

grid = RectilinearGrid(CPU(), size=(256, 32), extent=(500meters, 100meters), topology=(Bounded, Flat, Bounded))

biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid, open_bottom=true)

model = NonhydrostaticModel(; grid, biogeochemistry,
                              buoyancy=BuoyancyTracer(), tracers=:b,
                              advection=WENO(; grid),
                              closure = AnisotropicMinimumDissipation())

bᵢ(x, y, z) = ifelse(x < 250, 1e-4, 1e-3)

set!(model, b = bᵢ, N = 5.0, P = 0.1, Z = 0.1, T = 18.0)

simulation = Simulation(model; Δt = 1.0, stop_time = 3hours)

wizard = TimeStepWizard(cfl = 0.3, max_change = 1.5)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "buoyancy_front.jld2",
                                                       schedule = TimeInterval(1minute),
                                                       overwrite_existing = true)

run!(simulation)
```

We can then load the saved output and visualize it:

```@example quickstart
b = FieldTimeSeries("buoyancy_front.jld2", "b")
P = FieldTimeSeries("buoyancy_front.jld2", "P")

xb, yb, zb = nodes(b)
xP, yP, zP = nodes(P)

times = b.times

using CairoMakie

n = Observable(1)

b_lims = (minimum(b), maximum(b))
P_lims = (minimum(P), maximum(P))

bₙ = @lift interior(b[$n], :, 1, :)
Pₙ = @lift interior(P[$n], :, 1, :)

fig = Figure(resolution = (1200, 480), fontsize = 20)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 970)
ax1 = Axis(fig[1, 1]; title = "Buoyancy perturbation (m / s)", axis_kwargs...)
ax2 = Axis(fig[2, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)

hm1 = heatmap!(ax1, xb, zb, bₙ, colorrange = b_lims, colormap = :batlow)
hm2 = heatmap!(ax2, xP, zP, Pₙ, colorrange = P_lims, colormap = Reverse(:bamako))

Colorbar(fig[1, 2], hm1)
Colorbar(fig[2, 2], hm2)

record(fig, "buoyancy_front.mp4", 1:length(times)) do i
    n[] = i
end

nothing #hide
```

![buoyancy_front](buoyancy_front.mp4)

In this example `OceanBioME` is providing the `biogeochemistry` and the remainder is taken care of by `Oceananigans`. For comprehensive documentation of the physics modelling see [Oceananigans' Documentation](https://clima.github.io/OceananigansDocumentation/stable/); for biogeochemistry and other features we provide read below.

## Places to find OceanBioME information

* This documentation, which provides
    * example scripts,
    * explanations of model implementation methods,
    * details of currently implemented models, and
    * a library documenting all user-facing objects and functions.

* [Discussions on the OceanBioME github](https://github.com/OceanBioME/OceanBioME.jl/discussions)
  
    If you've got a question or something to talk about, don't hesitate to [start a new discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions/new?)!

* [Issues](https://github.com/OceanBioME/OceanBioME.jl//issues) and [pull requests](https://github.com/OceanBioME/OceanBioME.jl/pulls) also contain lots of information about problems we've found, solutions we're trying to implement, and ideas for the future.

## Getting in touch

Whether you need help getting started with OceanBioME, found a bug, want OceanBioME to be more expanded, or just want to chat about our project, you can:

* [Start a discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions). 
* [Open an issue](https://github.com/OceanBioME/OceanBioME.jl/issues). Issues are best if you think the OceanBioME source code needs attention: a bug, a sign error, an important missing feature, or a typo in the documentation.

## Citing

If you use OceanBioME as part of your research, teaching, or other activities, we would be grateful if you could cite our work and mention OceanBioME by name, as well as citing and [acknowledging Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/#Citing) as without them this package would not be possible.

We do not currently have a citation for OceanBioME so please reach out if you wish to cite it, and we will expedite the process of [making it citable](https://joss.theoj.org/about).

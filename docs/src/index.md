# `Ocean` `Bio`geochemical `M`odelling `E`nvironment - OceanBioME

```@meta
CurrentModule = OceanBioME
DocTestSetup = quote
    using OceanBioME
end
```

OceanBioME.jl aims to be a flexible and friendly ocean biogeochemical modelling environment for simulating coupled physics-biogeochemical models across a range of ocean scenarios. It is built within [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) which resolves the transport of the biogeochemical tracers and individuals and provides an excellent framework for implementing advanced features such as complex boundary conditions (e.g. sediment models).

OceanBioME.jl is currently developed primarily through the generous and greatly appreciated support of the [Center for Climate Repair at Cambridge](https://www.climaterepair.cam.ac.uk/) and the [Gordon and Betty Moore Foundation](https://www.moore.org/). 

OceanBioME.jl currently provides a core of several biogeochemical models (NPZD and [LOBSTER](https://doi.org/10.1029/2004JC002588), a medium complexity model, and [PISCES](https://doi.org/10.5194/gmd-8-2465-2015) in an early stage of testing), air-sea gas exchange models to provide appropriate top boundary conditions, and example sediment models for the benthic boundary.

As it was originally conceived for the study of submesoscale systems such as kelp forests we have also implemented a framework for integrating the growth of biological/active Lagrangian particles which move around and can interact with the (Eulerian) tracer fields - consuming nutrients and carbon dioxide, and depositing waste. Currently, included is a simple sugar kelp growth model which can be used in a variety of dynamical scenarios such as free floating or fixed to the ground, with more macroalgae models expected to be ready soon.

## Installation:

First, [download and install Julia](https://julialang.org/downloads/)

From the Julia prompt (REPL), type:
```julia
julia> using Pkg
julia> Pkg.add("OceanBioME")
```

## Running your first model
As a simple example lets run a Nutrient-Phytoplankton-Zooplankton-Detritus (NPZD) model in a two-dimensional simulation of a buoyancy front:
```julia
using OceanBioME, Oceananigans
using Oceananigans.Units

grid = RectilinearGrid(CPU(), size=(256, 32), extent=(500, 100), topology=(Bounded, Flat, Bounded))

biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid, open_bottom=true)

model = NonhydrostaticModel(; grid, biogeochemistry, buoyancy=BuoyancyTracer(), tracers=:b, advection=WENO(; grid), closure = AnisotropicMinimumDissipation())

bᵢ(x, y, z) = ifelse(x < 250, 1e-4, 1e-3)

set!(model, b = bᵢ, N = 5.0, P = 0.1, Z = 0.1, T = 18.0)

simulation = Simulation(model; Δt=1.0, stop_time=3hours)

wizard = TimeStepWizard(cfl=0.3, max_change=1.5)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers, filename = "buoyancy_front.jld2", schedule = TimeInterval(1minute), overwrite_existing=true)

run!(simulation)
```

<details>
<summary>We can then visualise this:</summary>

```julia
using CairoMakie
b = FieldTimeSeries("buoyancy_front.jld2", "b")
P = FieldTimeSeries("buoyancy_front.jld2", "P")

n = Observable(1)

b_lims = (minimum(b), maximum(b))
P_lims = (minimum(P), maximum(P))

b_plt = @lift b[1:grid.Nx, 1, 1:grid.Nz, $n]
P_plt = @lift P[1:grid.Nx, 1, 1:grid.Nz, $n]

fig = Figure(resolution = (1600, 160 * 4))

supertitle = Label(fig[0, :], "t = 0.0")

ax1 = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "z (m)", title = "Buouyancy pertubation (m / s)", width = 1400)
ax2 = Axis(fig[2, 1], xlabel = "x (m)", ylabel = "z (m)", title = "Phytoplankton concentration (mmol N / m³)", width = 1400)

hm1 = heatmap!(ax1, xnodes(Center, grid)[1:grid.Nx], znodes(Center, grid)[1:grid.Nz], b_plt, colorrange = b_lims, colormap = :batlow, interpolate=true)
hm2 = heatmap!(ax2, xnodes(Center, grid)[1:grid.Nx], znodes(Center, grid)[1:grid.Nz], P_plt, colorrange = P_lims, colormap = Reverse(:bamako), interpolate=true)

Colorbar(fig[1, 2], hm1)
Colorbar(fig[2, 2], hm2)

record(fig, "buoyancy_front.gif", 1:length(b.times)) do i
    n[] = i
    msg = string("Plotting frame ", i, " of ", length(b.times))
    print(msg * " \r")
    supertitle.text = "t=$(prettytime(b.times[i]))"
end
```
</details>

![buoyancy_front](https://user-images.githubusercontent.com/26657828/226373754-42c5c9ed-d7fc-450a-8346-a497a40fe0e2.gif)

In this example `OceanBioME` is providing the `biogeochemistry` and the remainder is taken care of by `Oceanaigans`. For comprehensive documentation of the physics modelling see [Oceanaigans' Documentation](https://clima.github.io/OceananigansDocumentation/stable/), and for biogeochemistry and other features we provide read below.

## Places to find OceanBioME information

* This documentation, which provides
    * example scripts,
    * explanations of model implementation methods,
    * details of currently implimented models, and
    * a library documenting all user-facing objects and functions.
* [Discussions on the OceanBioME github](https://github.com/OceanBioME/OceanBioME.jl/discussions)
  
    If you've got a question or something to talk about, don't hesitate to [start a new discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions/new?)!

* [Issues](https://github.com/OceanBioME/OceanBioME.jl//issues) and [pull requests](https://github.com/OceanBioME/OceanBioME.jl/pulls) also contain lots of information about problems we've found, solutions we're trying to implement, and ideas for the future.

## Getting in touch

Whether you need help getting started with OceanBioME, found a bug, want OceanBioME to be more expanded, or just want to chat about our project, you can:

* [Start a discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions). 
* [Open an issue](https://github.com/OceanBioME/OceanBioME.jl/issues). Issues are best if you think the OceanBioME source code needs attention: a bug, a sign error, an important missing feature, or a typo in this documentation.

## Citing

If you use OceanBioME as part of your research, teaching, or other activities, we would be grateful if you could cite our work and mention OceanBioME by name, as well as citing and [acknowledging Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/#Citing) as without them this package would not be possible.

We do not currently have a citation for OceanBioME so please reach out if you wish to cite it, and we will expedite the process of [making it citable](https://joss.theoj.org/about).

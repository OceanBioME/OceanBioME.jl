![](OceanBioME_headerbar.jpg?raw=true)
[![Documentation](https://img.shields.io/badge/documentation-stable%20release-blue?style=flat-square)](https://oceanbiome.github.io/OceanBioME.jl/stable/)
[![Documentation](https://img.shields.io/badge/documentation-dev%20release-orange?style=flat-square)](https://oceanbiome.github.io/OceanBioME.jl/dev/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://mit-license.org)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Ask us anything: discussion](https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg?style=flat-square)](https://github.com/OceanBioME/OceanBioME.jl/discussions)
[![GitHub tag (latest SemVer pre-release)](https://img.shields.io/github/v/tag/OceanBioME/OceanBioME.jl?include_prereleases&label=latest%20version&logo=github&sort=semver&style=flat-square)](https://github.com/OceanBioME/OceanBioME.jl/releases)

[![Testing](https://github.com/OceanBioME/OceanBioME.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/OceanBioME/OceanBioME.jl/actions/workflows/tests.yml)
[![Documentation](https://github.com/OceanBioME/OceanBioME.jl/actions/workflows/documentation.yml/badge.svg)](https://github.com/OceanBioME/OceanBioME.jl/actions/workflows/documentation.yml)
# *Ocean* *Bio*geochemical *M*odelling *E*nvironment

## Description
OceanBioME was developed with generous support from the Cambridge Centre for Climate Repair [CCRC](https://www.climaterepair.cam.ac.uk) and the Gordon and Betty Moore Foundation as a tool to study the effectiveness and impacts of ocean carbon dioxide removal (CDR) strategies.

OceanBioME is a flexible modelling environment written in Julia for modelling the coupled interactions between ocean biogeochemistry, carbonate chemistry, and physics. OceanBioME can be run as a stand-alone box model, or coupled with [Oceananigans.jl](https://github.com/cliMA/oceananigans.jl/) to run as a 1D column model or with 2 and 3D physics. 

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

<details>
<summary>We can then visualise this:</summary>

```julia
using CairoMakie

b = FieldTimeSeries("buoyancy_front.jld2", "b")
P = FieldTimeSeries("buoyancy_front.jld2", "P")

xb, yb, zb = nodes(b)
xP, yP, zP = nodes(P)

times = b.times

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

hm1 = heatmap!(ax1, xb, zb, bₙ, colorrange = b_lims, colormap = :batlow, interpolate=true)
hm2 = heatmap!(ax2, xP, zP, Pₙ, colorrange = P_lims, colormap = Reverse(:bamako), interpolate=true)

Colorbar(fig[1, 2], hm1)
Colorbar(fig[2, 2], hm2)

record(fig, "buoyancy_front.gif", 1:length(times)) do i
    @info string("Plotting frame ", i, " of ", length(times))
    n[] = i
end
```
</details>

![buoyancy_front](https://github.com/OceanBioME/OceanBioME.jl/assets/7112768/84f7f712-5648-4293-be18-608a4a3413ba)

In this example `OceanBioME` is providing the `biogeochemistry` and the remainder is taken care of by `Oceananigans`. For comprehensive documentation of the physics modelling see [Oceananigans' Documentation](https://clima.github.io/OceananigansDocumentation/stable/), and for biogeochemistry and other features we provide read below.

## Documentation

See the [documentation](https://oceanbiome.github.io/OceanBioME.jl) for full description and examples.

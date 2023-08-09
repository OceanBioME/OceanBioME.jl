![](OceanBioME_headerbar.jpg?raw=true)
[![Documentation](https://img.shields.io/badge/documentation-stable%20release-blue?style=flat-square)](https://oceanbiome.github.io/OceanBioME.jl/stable/)
[![Documentation](https://img.shields.io/badge/documentation-dev%20release-orange?style=flat-square)](https://oceanbiome.github.io/OceanBioME.jl/dev/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg?style=flat-square)](https://mit-license.org)
[![status](https://joss.theoj.org/papers/f991f4a8f0fc5ab1aba3dd1ac51f20bd/status.svg)](https://joss.theoj.org/papers/f991f4a8f0fc5ab1aba3dd1ac51f20bd)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![Ask us anything: discussion](https://img.shields.io/badge/Ask%20us-anything-1abc9c.svg?style=flat-square)](https://github.com/OceanBioME/OceanBioME.jl/discussions)
[![GitHub tag (latest SemVer pre-release)](https://img.shields.io/github/v/tag/OceanBioME/OceanBioME.jl?include_prereleases&label=latest%20version&logo=github&sort=semver&style=flat-square)](https://github.com/OceanBioME/OceanBioME.jl/releases)

[![Testing](https://github.com/OceanBioME/OceanBioME.jl/actions/workflows/tests.yml/badge.svg)](https://github.com/OceanBioME/OceanBioME.jl/actions/workflows/tests.yml)
[![codecov](https://codecov.io/gh/OceanBioME/OceanBioME.jl/branch/main/graph/badge.svg?token=3DIW4R7N3R)](https://codecov.io/gh/OceanBioME/OceanBioME.jl)
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
As a simple example lets run a Nutrient-Phytoplankton-Zooplankton-Detritus (NPZD) model in a two-dimensional simulation of a buoyancy front. This example requires Oceananigans, so we install that first:

```julia
using Pkg; Pkg.add("Oceananigans")

using OceanBioME, Oceananigans
using Oceananigans.Units

grid = RectilinearGrid(CPU(), size = (160, 32), extent = (10000meters, 500meters), topology = (Bounded, Flat, Bounded))

biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid) 

model = NonhydrostaticModel(; grid, biogeochemistry,
                              advection = WENO(; grid),
			                  closure = AnisotropicMinimumDissipation(),
			                  buoyancy = SeawaterBuoyancy(constant_salinity = true))

@inline front(x, z, μ, δ) = μ + δ * tanh((x - 7000 + 4 * z) / 500)

Pᵢ(x, y, z) = ifelse(z > -50, 0.03, 0.01) 
Nᵢ(x, y, z) = front(x, z, 2.5, -2)
Tᵢ(x, y, z) = front(x, z, 9, 0.05)

set!(model, N = Nᵢ, P = Pᵢ, Z = Pᵢ, T = Tᵢ)

simulation = Simulation(model; Δt = 50, stop_time = 4days)

simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "buoyancy_front.jld2",
                                                       schedule = TimeInterval(24minute),
                                                       overwrite_existing = true)

run!(simulation)
```

<details>
<summary>We can then visualise this:</summary>

```julia
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

fig = Figure(resolution = (1000, 520), fontsize = 20)

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

record(fig, "buoyancy_front.gif", 1:length(times)) do i
    n[] = i
end
```
</details>

![buoyancy_front](https://github.com/OceanBioME/OceanBioME.jl/assets/26657828/b5cb0bba-b3a5-4600-a7bb-bccb825024a2)

In this example `OceanBioME` is providing the `biogeochemistry` and the remainder is taken care of by `Oceananigans`.
For comprehensive documentation of the physics modelling see
[Oceananigans' Documentation](https://clima.github.io/OceananigansDocumentation/stable/), and for
biogeochemistry and other features we provide read below.

## Using GPU

To run the same example on a GPU we just need to construct the `grid` on the GPU; the rest is taken care of!

Just replace `CPU()` with `GPU()` in the grid construction with everything else left unchanged:

```julia
grid = RectilinearGrid(GPU(), size = (256, 32), extent = (500meters, 100meters), topology = (Bounded, Flat, Bounded))
```

## Documentation

See the [documentation](https://oceanbiome.github.io/OceanBioME.jl) for full description of the software
package and more examples.

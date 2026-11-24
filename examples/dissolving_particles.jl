# # [Sinking and dissolving particles in two-dimensional turbulence](@id dissolving-particles-example)
#
# In this example we release mineral particles at the surface of a two-dimensional turbulent flow,
# as might be done for ocean alkalinity enhancement. The particles sink and dissolve, releasing
# alkalinity into the water. Large particles sink out of the eddies before they dissolve, while
# small particles dissolve near the surface.

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, CairoMakie, JLD2"
# ```

# ## Model setup
using OceanBioME, Oceananigans, Oceananigans.Units
using CairoMakie, JLD2, Printf, Random

Random.seed!(1234)
nothing #hide

# We use a two-dimensional domain which is 64 m wide and deep, and periodic in x
grid = RectilinearGrid(size = (128, 128), x = (0, 64), z = (-64, 0),
                       topology = (Periodic, Flat, Bounded))

# ## Particles
# [`DissolvingParticles`](@ref) are spheres whose radius shrinks at a `dissolution_rate`, and which
# sink at a `sinking_speed`. By default the radius shrinks at a constant rate and the particles sink
# according to Stokes' law, so they slow down as they dissolve. Both can instead be functions of the
# particle properties and the model fields at the particle, see [`DissolvingMineral`](@ref).
n = 1000

particles = DissolvingParticles(n; grid)

# The particles start with zero radius, so they do nothing until we release them.
#
# We couple the particles to a carbonate system without any biology, which has the dissolved
# inorganic carbon (`DIC`) and alkalinity (`Alk`) tracers that the particles release
biogeochemistry = NutrientsPlanktonDetritus(grid; inorganic_carbon = CarbonateSystem(), particles)

model = NonhydrostaticModel(grid; biogeochemistry, advection = WENO(order = 5))

# ## Spin up
# We start from random velocities and let the turbulence develop for 30 minutes
Alk₀ = 2300 # mmol/m³

uᵢ(x, z) = 0.02 * randn() # m/s

set!(model, u = uᵢ, w = uᵢ, DIC = 2000, Alk = Alk₀)

simulation = Simulation(model, Δt = 1, stop_time = 30minutes)

conjure_time_step_wizard!(simulation, cfl = 0.5, max_Δt = 10)

progress(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                        iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(sim.run_wall_time))

add_callback!(simulation, progress, TimeInterval(30minutes))

run!(simulation)

# ## Release
# Now we release the particles in a 16 m wide patch at the surface, with radii between 20 and 80 μm.
#
# Each computational particle stands for many real particles, the number of which is set by the
# `scalefactors`. Here we release 1 kg of mineral per metre in the (flat) y direction.
r₀ = 20e-6 .+ 60e-6 * rand(n)
x₀ = 24 .+ 16 * rand(n)
z₀ = - rand(n)

mineral = particles.biogeochemistry

scalefactor = 1 / (mineral.density * 4π / 3 * sum(r₀.^3)) # 1 kg / total particle mass

set!(particles, r = r₀, x = x₀, z = z₀, scalefactors = scalefactor)

# We save the alkalinity, its horizontal average and integral, and the particle positions and radii
Alk = model.tracers.Alk

filename = "dissolving_particles"

simulation.output_writers[:alkalinity] = JLD2Writer(model, (; Alk, mean_Alk = Average(Alk, dims = 1), total_Alk = Integral(Alk));
                                                    filename = "$filename.jld2",
                                                    schedule = TimeInterval(1minute),
                                                    overwrite_files = true)

particle_state(model) = (; x = Array(particles.x), z = Array(particles.z), r = Array(particles.fields.r))

simulation.output_writers[:particles] = JLD2Writer(model, (; particles = particle_state);
                                                   filename = "$(filename)_particles.jld2",
                                                   schedule = TimeInterval(1minute),
                                                   overwrite_files = true)

# and run for long enough that the largest particles dissolve
simulation.stop_time = 3hours

run!(simulation)

# ## Results
# First we load the alkalinity and particle output
Alk = FieldTimeSeries("$filename.jld2", "Alk")
mean_Alk = FieldTimeSeries("$filename.jld2", "mean_Alk")
total_Alk = FieldTimeSeries("$filename.jld2", "total_Alk")

times = Alk.times

file = jldopen("$(filename)_particles.jld2")

states = [file["timeseries/particles/$iteration"] for iteration in keys(file["timeseries/t"])]

close(file)

x, z, r = (stack(getproperty.(states, name)) for name in (:x, :z, :r))
nothing #hide

# and animate the particles sinking and dissolving through the turbulence, with each particle drawn
# at its current size (magnified)
frame = Observable(1)

Alkₙ = @lift Alk[$frame] - Alk₀
positions = @lift Point2f.(x[:, $frame], z[:, $frame])
sizes = @lift 2e5 * r[:, $frame]
title = @lift prettytime(times[$frame] - times[1]) * " after release"

fig = Figure(size = (700, 600), fontsize = 20)

ax = Axis(fig[1, 1]; title, xlabel = "x (m)", ylabel = "z (m)", aspect = DataAspect())

hm = heatmap!(ax, Alkₙ, colormap = Reverse(:deep), colorrange = (0, 80))
scatter!(ax, positions, markersize = sizes, color = :white)

Colorbar(fig[1, 2], hm, label = "Alkalinity increase (mmol/m³)")

record(fig, "$filename.mp4", eachindex(times), framerate = 12) do i
    frame[] = i
end
nothing #hide

# ![](dissolving_particles.mp4)

# To see how the turbulence changes where the particles dissolve, we compare the depth at which each
# particle finished dissolving with the depth it would reach in still water. A particle sinking at
# the Stokes' law speed, ``w ∝ r²``, while its radius shrinks at a constant rate ``k``, sinks
# ``|w(r₀)| r₀ / 3k`` before it has dissolved.
dissolution_depth = [-z[p, findfirst(iszero, r[p, :])] for p in 1:n]

k = mineral.dissolution_rate

still_water_depth(r) = - mineral.sinking_speed(mineral, (; r), nothing) * r / 3k
nothing #hide

# We also check that the alkalinity lost from the particles is gained by the water
θ = mineral.released_tracers.Alk # mol Alk / mol mineral

mineral_concentration = 1000 * mineral.density / mineral.molar_mass # mmol / m³

particle_alkalinity = scalefactor * θ * mineral_concentration * 4π / 3 * vec(sum(r.^3, dims = 1))
dissolved_alkalinity = interior(total_Alk, 1, 1, 1, :) .- total_Alk[1, 1, 1, 1]
nothing #hide

# and plot these along with the particle radii and the horizontally averaged alkalinity increase
hours_after_release = (times .- times[1]) / hour

fig = Figure(size = (1200, 900), fontsize = 20)

ax = Axis(fig[1, 1], xlabel = "Initial radius (μm)", ylabel = "Depth (m)", title = "Depth when fully dissolved",
          yreversed = true, limits = (nothing, (0, 64)))

scatter!(ax, 1e6 * r₀, dissolution_depth, color = (:black, 0.3), label = "Turbulence")
lines!(ax, 20:80, still_water_depth.(1e-6 * (20:80)), color = :red, linewidth = 3, label = "Still water")
axislegend(ax, position = :lb)

ax = Axis(fig[1, 3], xlabel = "Time (hours)", ylabel = "Radius (μm)", title = "Particle radius")

for p in 1:10:n
    lines!(ax, hours_after_release, 1e6 * r[p, :], color = r₀[p], colormap = :viridis, colorrange = (20e-6, 80e-6))
end

ax = Axis(fig[2, 1], xlabel = "Time (hours)", ylabel = "z (m)", title = "Horizontally averaged alkalinity increase")

hm = heatmap!(ax, hours_after_release, znodes(mean_Alk), interior(mean_Alk, 1, 1, :, :)' .- Alk₀, colormap = Reverse(:deep))
Colorbar(fig[2, 2], hm, label = "mmol/m³")

ax = Axis(fig[2, 3], xlabel = "Time (hours)", ylabel = "Fraction of alkalinity", title = "Alkalinity budget")

lines!(ax, hours_after_release, particle_alkalinity / particle_alkalinity[1], linewidth = 3, label = "In particles")
lines!(ax, hours_after_release, dissolved_alkalinity / particle_alkalinity[1], linewidth = 3, label = "Dissolved")
lines!(ax, hours_after_release, (particle_alkalinity .+ dissolved_alkalinity) / particle_alkalinity[1],
       linewidth = 3, linestyle = :dash, color = :black, label = "Total")
axislegend(ax, position = :rc)

fig

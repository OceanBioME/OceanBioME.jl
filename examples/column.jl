# # One dimensional column example
# In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution. This demonstrates:
# - How to setup OceanBioME's biogeochemical models
# - How to visualise results
# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, Printf, CairoMakie"
# ```

# ## Model setup
# We load the packages and choose the default LOBSTER parameter set
using OceanBioME, Oceananigans, Printf
using OceanBioME.SLatissimaModel: SLatissima
using Oceananigans.Units

const year = years = 365days

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 +exp(-(t - 100days) / (5days)))) * (1 / (1 + exp((t - 330days) / (25days))))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(year-eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t))/10)) / 2 + 1e-4

@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50day) + 10

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
Lx, Ly = 20meters, 20meters
grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200meters)) 

# Specify the boundary conditions for DIC and O₂ based on the air-sea CO₂ and O₂ flux
CO₂_flux = GasExchange(; gas = :CO₂, temperature = temp, salinity = (args...) -> 35)

model = NonhydrostaticModel(; grid,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                              biogeochemistry = LOBSTER(; grid,
                                                          surface_phytosynthetically_active_radiation = PAR⁰,
                                                          carbonates = true,
                                                          advection_schemes = (sPOM = WENO(grid), bPOM = WENO(grid))),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ),
                              advection = nothing)

set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt=3minutes, stop_time=100days) 

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))      
                                                                  
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))

filename = "column"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)


#simulation.output_writers[:particles] = JLD2OutputWriter(model, (; particles), filename = "$(filename)_particles.jld2", schedule = TimeInterval(1day), overwrite_existing = true)

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())

# ## Run!
# Finally we run the simulation
run!(simulation)

# Now we can visulise the results with some post processing to diagnose the air-sea CO₂ flux - hopefully this looks different to the example without kelp!

   P = FieldTimeSeries("$filename.jld2", "P")
 NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
   Z = FieldTimeSeries("$filename.jld2", "Z")
sPOM = FieldTimeSeries("$filename.jld2", "sPOM") 
bPOM = FieldTimeSeries("$filename.jld2", "bPOM")
 DIC = FieldTimeSeries("$filename.jld2", "DIC")
 Alk = FieldTimeSeries("$filename.jld2", "Alk")

x, y, z = nodes(P)
times = P.times

air_sea_CO₂_flux = zeros(length(times))
carbon_export = zeros(length(times))
for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = CO₂_flux.condition.parameters(0.0, 0.0, t, DIC[1, 1, 50, i], Alk[1, 1, 50, i], temp(1, 1, 0, t), 35)
    carbon_export[i] = (sPOM[1, 1, end-20, i] * model.biogeochemistry.sinking_velocities.sPOM.w[1, 1, end-20] + bPOM[1, 1, end-20, i] * model.biogeochemistry.sinking_velocities.bPOM.w[1, 1, end-20]) * model.biogeochemistry.organic_redfield
end

using CairoMakie

fig = Figure(resolution = (1920, 1050), fontsize=24)

axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0, times[end] / days), (-150, 0)))
heatmap_kwargs = (interpolate = true, colormap = :batlow)

axP = Axis(fig[1, 1:2]; title="Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hmP = heatmap!(times / days, z, interior(P, 1, 1, :, :)'; heatmap_kwargs...)
cbP = Colorbar(fig[1, 3], hmP)

axNO₃ = Axis(fig[1, 4:5]; title="Nitrate concentration (mmol N / m³)", axis_kwargs...)
hmNO₃ = heatmap!(times / days, z, interior(NO₃, 1, 1, :, :)'; heatmap_kwargs...)
cbNO₃ = Colorbar(fig[1, 6], hmNO₃)

axZ = Axis(fig[2, 1:2]; title="Zooplankton concentration (mmol N / m³)", axis_kwargs...)
hmZ = heatmap!(times / days, z, interior(Z, 1, 1, :, :)'; heatmap_kwargs...)
cbZ = Colorbar(fig[2, 3], hmZ)

axD = Axis(fig[2, 4:5]; title="Detritus concentration (mmol N / m³)", axis_kwargs...)
hmD = heatmap!(times / days, z, interior(sPOM, 1, 1, :, :)' .+ interior(bPOM, 1, 1, :, :)'; heatmap_kwargs...)
cbD = Colorbar(fig[2, 6], hmD)

axfDIC = Axis(fig[3, 1:4], xlabel = "Time (days)", ylabel = "Flux (kgCO₂ / m² / year)", title = "Air-sea CO₂ flux and Sinking")
hmfDIC = lines!(times / days, air_sea_CO₂_flux * (12 + 16 * 2) * year / 1e6, linewidth=3, label = "Air-sea flux")
hmfExp = lines!(times / days, carbon_export    * (12 + 16 * 2) * year / 1e6, linewidth=3, label = "Sinking export")

fig[3, 5:6] = Legend(fig, axfDIC, "", framevisible = false)

current_figure() # hide
fig

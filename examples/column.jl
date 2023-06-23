# # One dimensional column example
# In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution. This demonstrates:
# - How to setup OceanBioME's biogeochemical models
# - How to visualise results
# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, CairoMakie"
# ```

# ## Model setup
# We load the packages and choose the default LOBSTER parameter set
using OceanBioME, Oceananigans, Printf
using OceanBioME.SLatissimaModel: SLatissima
using Oceananigans.Units

const year = years = 365days
nothing #hide

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4

@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50days) + 10
nothing #hide

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
                                                          carbonates = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ),
                              advection = nothing)

set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)

simulation = Simulation(model, Δt = 3minutes, stop_time = 100days)

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))      
                                                                  
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))

filename = "column"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields),
                                                        filename = "$filename.jld2",
                                                        schedule = TimeInterval(1day),
                                                        overwrite_existing = true)

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())
nothing #hide

# ## Run!
# We are ready to run the simulation
run!(simulation)

# ## Load saved output
# Now we can load the results and do some post processing to diagnose the air-sea CO₂ flux. Hopefully, this looks different to the example without kelp!

   P = FieldTimeSeries("$filename.jld2", "P")
 NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
   Z = FieldTimeSeries("$filename.jld2", "Z")
sPOM = FieldTimeSeries("$filename.jld2", "sPOM")
bPOM = FieldTimeSeries("$filename.jld2", "bPOM")
 DIC = FieldTimeSeries("$filename.jld2", "DIC")
 Alk = FieldTimeSeries("$filename.jld2", "Alk")

x, y, z = nodes(P)
times = P.times

# We compute the  air-sea CO₂ flux at the surface (corresponding to vertical index `k = grid.Nz`) and
# the carbon export by computing how much carbon sinks below some arbirtrary depth; here we use depth 
# that corresponds to `k = grid.Nz - 20`.
air_sea_CO₂_flux = zeros(length(times))
carbon_export = zeros(length(times))

for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = CO₂_flux.condition.parameters(0.0, 0.0, t, DIC[1, 1, grid.Nz, i], Alk[1, 1, grid.Nz, i], temp(1, 1, 0, t), 35)
    carbon_export[i] = (sPOM[1, 1, grid.Nz-20, i] * model.biogeochemistry.sinking_velocities.sPOM.w[1, 1, grid.Nz-20] +
                        bPOM[1, 1, grid.Nz-20, i] * model.biogeochemistry.sinking_velocities.bPOM.w[1, 1, grid.Nz-20]) * model.biogeochemistry.organic_redfield
end

# Both `air_sea_CO₂_flux` and `carbon_export` are in units `mmol CO₂ / (m² s)`.

# ## Plot
# Finally, we plot!

using CairoMakie

fig = Figure(resolution = (1000, 1500), fontsize = 20)

axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0, times[end] / days), (-150meters, 0)))

axP = Axis(fig[1, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hmP = heatmap!(times / days, z, interior(P, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[1, 2], hmP)

axNO₃ = Axis(fig[2, 1]; title = "Nitrate concentration (mmol N / m³)", axis_kwargs...)
hmNO₃ = heatmap!(times / days, z, interior(NO₃, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[2, 2], hmNO₃)

axZ = Axis(fig[3, 1]; title = "Zooplankton concentration (mmol N / m³)", axis_kwargs...)
hmZ = heatmap!(times / days, z, interior(Z, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[3, 2], hmZ)

axD = Axis(fig[4, 1]; title = "Detritus concentration (mmol N / m³)", axis_kwargs...)
hmD = heatmap!(times / days, z, interior(sPOM, 1, 1, :, :)' .+ interior(bPOM, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[4, 2], hmD)

CO₂_molar_mass = (12 + 2 * 16) * 1e-3 # kg / mol

axfDIC = Axis(fig[5, 1], xlabel = "Time (days)", ylabel = "Flux (kgCO₂/m²/year)",
                         title = "Air-sea CO₂ flux and Sinking", limits = ((0, times[end] / days), nothing))
lines!(axfDIC, times / days, air_sea_CO₂_flux / 1e3 * CO₂_molar_mass * year, linewidth = 3, label = "Air-sea flux")
lines!(axfDIC, times / days, carbon_export / 1e3    * CO₂_molar_mass * year, linewidth = 3, label = "Sinking export")
Legend(fig[5, 2], axfDIC, framevisible = false)

fig

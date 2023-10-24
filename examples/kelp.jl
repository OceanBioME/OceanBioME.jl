# # Simple active particle example
# Here, we setup a simple 1D column example with the [LOBSTER](@ref LOBSTER) biogeochemical model and active particles modelling the growth of sugar kelp.
# This example demonstrates:
# - How to setup OceanBioME's biogeochemical models
# - How to add biologically active particles which interact with the biodeochemical model
# - How to visualise results

# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first.

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg "add OceanBioME, Oceananigans, CairoMakie, JLD2"
# ```

# ## Model setup
# First load the required packages
using OceanBioME, Oceananigans, Printf
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units
using Oceananigans.Architectures: arch_array

const year = years = 365days # just for these idealised cases
nothing #hide

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic).

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days) ^ 2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t))/10)) / 2 + 1e-4

@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50day) + 10
nothing #hide

architecture = CPU()

# ## Grid and PAR field
# Define the grid and an extra Oceananigans' field that the PAR will be stored in
Lx, Ly = 20meters, 20meters
grid = RectilinearGrid(architecture, size=(1, 1, 50), extent=(Lx, Ly, 200)) 

# Specify the boundary conditions for DIC and O₂ based on the air-sea CO₂ and O₂ flux
CO₂_flux = GasExchange(; gas = :CO₂)

clock = Clock(; time = 0.0)
T = FunctionField{Center, Center, Center}(temp, grid; clock)
S = ConstantField(35)

# ## Kelp Particle setup
@info "Setting up kelp particles"

n = 5 # number of kelp fronds
z₀ = [-21:5:-1;] * 1.0 # depth of kelp fronds

particles = SLatissima(; architecture, 
                         x = arch_array(architecture, ones(n) * Lx / 2), 
                         y = arch_array(architecture, ones(n) * Ly / 2), 
                         z = arch_array(architecture, z₀), 
                         A = arch_array(architecture, ones(n) * 10.0),
                         latitude = 57.5,
                         scalefactor = 500.0)

# ## Setup BGC model
biogeochemistry = LOBSTER(; grid,
                            surface_phytosynthetically_active_radiation = PAR⁰,
                            carbonates = true,
                            variable_redfield = true,
                            scale_negatives = true,
                            particles)

model = NonhydrostaticModel(; grid,
                              clock,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                              biogeochemistry,
                              auxiliary_fields = (; T, S))

set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implementation page)

simulation = Simulation(model, Δt = 3minutes, stop_time = 100days) 

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))      
                                                                  
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))

filename = "kelp"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, model.tracers,
                                                        filename = "$filename.jld2",
                                                        schedule = TimeInterval(1day),
                                                        overwrite_existing = true)

simulation.output_writers[:particles] = JLD2OutputWriter(model, (; particles),
                                                         filename = "$(filename)_particles.jld2",
                                                         schedule = TimeInterval(1day),
                                                         overwrite_existing = true)

nothing #hide

# ## Run!
# Finally we run the simulation
run!(simulation)

# Now we can visualise the results with some post processing to diagnose the air-sea CO₂ flux - hopefully this looks different to the example without kelp!

   P = FieldTimeSeries("$filename.jld2", "P")
 NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
   Z = FieldTimeSeries("$filename.jld2", "Z")
sPOC = FieldTimeSeries("$filename.jld2", "sPOC") 
bPOC = FieldTimeSeries("$filename.jld2", "bPOC")
 DIC = FieldTimeSeries("$filename.jld2", "DIC")
 Alk = FieldTimeSeries("$filename.jld2", "Alk")

x, y, z = nodes(P)
times = P.times
nothing #hide

# We compute the  air-sea CO₂ flux at the surface (corresponding to vertical index `k = grid.Nz`) and
# the carbon export by computing how much carbon sinks below some arbirtrary depth; here we use depth 
# that corresponds to `k = grid.Nz - 20`.
air_sea_CO₂_flux = zeros(length(times))
carbon_export = zeros(length(times))

using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity

for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = CO₂_flux.condition.func(0.0, 0.0, t, DIC[1, 1, grid.Nz, i], Alk[1, 1, grid.Nz, i], temp(1, 1, 0, t), 35)
    carbon_export[i] = sPOC[1, 1, grid.Nz-20, i] * biogeochemical_drift_velocity(biogeochemistry, Val(:sPOC)).w[1, 1, grid.Nz-20] +
                       bPOC[1, 1, grid.Nz-20, i] * biogeochemical_drift_velocity(biogeochemistry, Val(:bPOC)).w[1, 1, grid.Nz-20]
end

# Both `air_sea_CO₂_flux` and `carbon_export` are in units `mmol CO₂ / (m² s)`.

using CairoMakie

fig = Figure(resolution = (1000, 1500), fontsize = 20)

axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0, times[end] / days), (-85meters, 0)))

axP = Axis(fig[1, 1]; title = "Phytoplankton concentration (mmol N/m³)", axis_kwargs...)
hmP = heatmap!(times / days, z, interior(P, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[1, 2], hmP)

axNO₃ = Axis(fig[2, 1]; title = "Nitrate concentration (mmol N/m³)", axis_kwargs...)
hmNO₃ = heatmap!(times / days, z, interior(NO₃, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[2, 2], hmNO₃)

axZ = Axis(fig[3, 1]; title = "Zooplankton concentration (mmol N/m³)", axis_kwargs...)
hmZ = heatmap!(times / days, z, interior(Z, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[3, 2], hmZ)

axD = Axis(fig[4, 1]; title = "Detritus concentration (mmol C/m³)", axis_kwargs...)
hmD = heatmap!(times / days, z, interior(sPOC, 1, 1, :, :)' .+ interior(bPOC, 1, 1, :, :)', colormap = :batlow)
Colorbar(fig[4, 2], hmD)

CO₂_molar_mass = (12 + 2 * 16) * 1e-3 # kg / mol

axfDIC = Axis(fig[5, 1], xlabel = "Time (days)", ylabel = "Flux (kgCO₂/m²/year)",
                         title = "Air-sea CO₂ flux and Sinking", limits = ((0, times[end] / days), nothing))
lines!(axfDIC, times / days, air_sea_CO₂_flux / 1e3 * CO₂_molar_mass * year, linewidth = 3, label = "Air-sea flux")
lines!(axfDIC, times / days, carbon_export / 1e3    * CO₂_molar_mass * year, linewidth = 3, label = "Sinking export")
Legend(fig[5, 2], axfDIC, framevisible = false)

fig

# We can also have a look at how the kelp particles evolve
using JLD2

file = jldopen("$(filename)_particles.jld2")

iterations = keys(file["timeseries/t"])

A, N, C = ntuple(n -> zeros(5, length(iterations)), 3)

times = zeros(length(iterations))

for (i, iter) in enumerate(iterations)
    particles = file["timeseries/particles/$iter"]
    A[:, i] = particles.A
    N[:, i] = particles.N
    C[:, i] = particles.C

    times[i] = file["timeseries/t/$iter"]
end

fig = Figure(resolution = (1000, 800), fontsize = 20)

axis_kwargs = (xlabel = "Time (days)", limits = ((0, times[end] / days), nothing))

ax1 = Axis(fig[1, 1]; ylabel = "Frond area (dm²)", axis_kwargs...)
[lines!(ax1, times / day, A[n, :], linewidth = 3) for n in 1:5]

ax2 = Axis(fig[2, 1]; ylabel = "Total Carbon (gC)", axis_kwargs...)
[lines!(ax2, times / day, (@. A * (N + particles.structural_nitrogen) * particles.structural_dry_weight_per_area)[n, :], linewidth = 3) for n in 1:5]

ax3 = Axis(fig[3, 1]; ylabel = "Total Nitrogen (gN)", axis_kwargs...)
[lines!(ax3, times / day, (@. A * (C + particles.structural_carbon) * particles.structural_dry_weight_per_area)[n, :], linewidth = 3) for n in 1:5]

fig

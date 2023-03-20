# # Simple active particle example
# In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and active particles modelling the growth of sugar kelp. This demonstrates:
# - How to setup OceanBioME's biogeochemical models
# - How to add biologically active particles which interact with the biodeochemical model
# - How to visualise results

# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg "add OceanBioME, Oceananigans, Printf, CairoMakie"
# ```

# ## Model setup
# First load the required packages
using OceanBioME, Oceananigans, Printf
using Oceananigans.Units

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic)

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / 365days))*(1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days) ^ 2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, 365days) * (1 / (1 +exp(-(t - 100days) / (5days)))) * (1 / (1 + exp((t - 330days) / (25days))))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(365days-eps(365days)) * exp(-mod(t, 365days) / 25days) - fmld1(mod(t, 365days))))

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t))/10)) / 2 + 1e-4

@inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50day) + 10

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
Lx, Ly = 20, 20
grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200)) 

# Specify the boundary conditions for DIC and O₂ based on the air-sea CO₂ and O₂ flux
CO₂_flux = GasExchange(; gas = :CO₂, temperature = temp, salinity = (args...) -> 35)

# ## Kelp Particle setup
@info "Setting up kelp particles"
n = 5 # number of kelp fronds
z₀ = [-21:5:-1;] * 1.0 # depth of kelp fronds

particles = SLatissima(; x = ones(n) * Lx / 2, y = ones(n) * Ly / 2, z = z₀, 
                         A = ones(n) * 10.0,
                         latitude = 57.5,
                         scalefactor = 500.0,
                         pescribed_temperature = temp)

# ## Setup BGC model
biogeochemistry = LOBSTER(; grid,
                            surface_phytosynthetically_active_radiation = PAR⁰,
                            carbonates = true,
                            variable_redfield = true)

model = NonhydrostaticModel(; grid,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                              biogeochemistry,
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ),
                              advection = nothing)


set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Show the progress of the simulation
# - Store the model and particles output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implementation page)

simulation = Simulation(model, Δt = 10minutes, stop_time=100days) 

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))      
                                                                  
simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))

filename = "kelp"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)

simulation.output_writers[:particles] = JLD2OutputWriter(model, (; particles), filename = "$(filename)_particles.jld2", schedule = TimeInterval(1day), overwrite_existing = true)

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())

wizard = TimeStepWizard(cfl = 0.1, diffusive_cfl = 0.1, max_change = 2.0, min_change = 0.5, cell_diffusion_timescale = column_diffusion_timescale, cell_advection_timescale = column_advection_timescale)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

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

air_sea_CO₂_flux = zeros(size(P)[4])
carbon_export = zeros(size(P)[4])
for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = CO₂_flux.condition.parameters(0.0, 0.0, t, DIC[1, 1, 50, i], Alk[1, 1, 50, i], temp(1, 1, 0, t), 35)
    carbon_export[i] = sPOC[1, 1, end - 20, i] * model.biogeochemistry.sinking_velocities.sPOM.w[1, 1, end - 20] + bPOC[1, 1, end - 20, i] * model.biogeochemistry.sinking_velocities.bPOM.w[1, 1, end - 20]
end

using CairoMakie
f = Figure(resolution = (1920, 1050))

axP = Axis(f[1, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Phytoplankton concentration (mmol N/m³)")
hmP = heatmap!(times./days, float.(z[end-23:end]), float.(P[1, 1, end-23:end, 1:end])', interpolate=true, colormap=:batlow)
cbP = Colorbar(f[1, 3], hmP)

axNO₃ = Axis(f[1, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Nitrate concentration (mmol N/m³)")
hmNO₃ = heatmap!(times./days, float.(z[end-23:end]), float.(NO₃[1, 1, end-23:end, 1:end])', interpolate=true, colormap=:batlow)
cbNO₃ = Colorbar(f[1, 6], hmNO₃)

axZ = Axis(f[2, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Zooplankton concentration (mmol N/m³)")
hmZ = heatmap!(times./days, float.(z[end-23:end]), float.(Z[1, 1, end-23:end, 1:end])', interpolate=true, colormap=:batlow)
cbZ = Colorbar(f[2, 3], hmZ)

axD = Axis(f[2, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Detritus concentration (mmol C/m³)")
hmD = heatmap!(times./days, float.(z[end-23:end]), float.(sPOC[1, 1, end-23:end, 1:end])' .+ float.(bPOC[1, 1, end-23:end, 1:end])', interpolate=true, colormap=:batlow)
cbD = Colorbar(f[2, 6], hmD)

axfDIC = Axis(f[3, 1:4], xlabel="Time (days)", title="Air-sea CO₂ flux and Sinking", ylabel="Flux (kgCO₂/m²/year)")
hmfDIC = lines!(times ./ days, air_sea_CO₂_flux .* (12 + 16 * 2) .* year /(1000 * 1000), label="Air-sea flux")
hmfExp = lines!(times ./ days, carbon_export .* (12 + 16 * 2) .* year / (1000 * 1000), label="Sinking export")

f[3, 5] = Legend(f, axfDIC, "", framevisible = false)

save("$(filename).png", f)

# ![Results-bgc](kelp.png)

# We can also have a look at how the kelp particles evolve
using JLD2

file = jldopen("$(filename)_particles.jld2")

iterations = keys(file["timeseries/t"])

A, N, C = ntuple(n -> ones(5, length(iterations)) .* NaN, 3)
times = ones(length(iterations)) .* NaN

for (i, iter) in enumerate(iterations)
    particles = file["timeseries/particles/$iter"]
    A[:, i] = particles.A
    N[:, i] = particles.N
    C[:, i] = particles.C

    times[i] = file["timeseries/t/$iter"]
end

fig = Figure(resolution = (1600, 500))

ax1 = Axis(fig[1, 1], ylabel = "Frond area (dm²)", xlabel = "Time (days)")
[lines!(ax1, times./day, A[n, :]) for n in 1:5]

ax2 = Axis(fig[1, 2], ylabel = "Total Carbon (gC)", xlabel = "Time (days)")
[lines!(ax2, times./day, (A .* (N .+ particles.structural_nitrogen) .* particles.structural_dry_weight_per_area)[n, :]) for n in 1:5]

ax3 = Axis(fig[1, 3], ylabel = "Total Nitrogen (gN)", xlabel = "Time (days)")
[lines!(ax3, times./day, (A .* (C .+ particles.structural_carbon) .* particles.structural_dry_weight_per_area)[n, :]) for n in 1:5]

save("kelp_particles.png", fig)
# ![Results](kelp_particles.png)=#

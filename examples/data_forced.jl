# # One dimensional column forced by external data with carbonate chemistry
# In this example we will setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution. This demonstraits:
# - How to setup OceanBioME's biogeochemical models
# - How to load external forcing data
# - How to run with optional tracer sets such as carbonate chemistry
# - How to setup a non-uniform grid for better near surface resolution
# - How to visulise results

# This is forced by mixing layer depth and surface photosynthetically available radiation (PAR) data from the Mercator Ocean model and NASA VIIRS observations

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, Printf, Plots, GLMakie, NetCDF, JLD2, DataDeps, Interpolations"
# ```

# ## Model setup
# We load the packages and choose the default LOBSTER parameter set

using Oceananigans, Random, Printf, NetCDF, Interpolations, DataDeps
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Operators: ∂zᶜᶜᶜ
using OceanBioME 

# ## Load external forcing data
# Loading the forcing data from our online copy
dd = DataDep(
    "example_data",
    "example data from subpolar re analysis and observational products", 
    "https://github.com/OceanBioME/OceanBioME_example_data/raw/main/subpolar.nc"
)
register(dd)
filename = datadep"example_data/subpolar.nc"
times = ncread(filename, "time")
temp = ncread(filename, "temp")
salinity = ncread(filename, "so")
mld = ncread(filename, "mld")
par = ncread(filename, "par")

temperature_itp = LinearInterpolation(times, temp) 
salinity_itp = LinearInterpolation(times, salinity) 
mld_itp = LinearInterpolation(times, mld) 
PAR_itp = LinearInterpolation(times, par)

t_function(x, y, z, t) = temperature_itp(mod(t, 364days))
s_function(x, y, z, t) = salinity_itp(mod(t, 364days))
surface_PAR(x, y, t) = PAR_itp(mod(t, 364days))
κₜ(x, y, z, t) = 2e-2*max(1-(z+mld_itp(mod(t,364days))/2)^2/(mld_itp(mod(t,364days))/2)^2,0)+1e-4

# ## Grid and PAR field
# Define the grid (in this case a non uniform grid for better resolution near the surface) and an extra Oceananigans field for the PAR to be stored in
Nz = 33
Lz = 600
refinement = 10 
stretching = 5.754
h(k) = (k - 1) / Nz
ζ₀(k) = 1 + (h(k) - 1) / refinement
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(size = (1, 1, Nz), x = (0, 20), y = (0, 20), z = z_faces)        
PAR = CenterField(grid)

# ## Biogeochemical and Oceananigans model
# Here we instantiate the LOBSTER model with carbonate chemistry and a surface flux of DIC (CO₂)
CO₂_flux = GasExchange(; gas = :CO₂, temperature = t_function, salinity = s_function)
model = NonhydrostaticModel(; grid,
                              closure = ScalarDiffusivity(ν=κₜ, κ=κₜ), 
                              biogeochemistry = LOBSTER(; grid,
                                                          surface_phytosynthetically_active_radiation = surface_PAR,
                                                          carbonates = true),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), ),
                              auxiliary_fields = (; PAR))

set!(model, P=0.03, Z=0.03, NO₃=11.0, NH₄=0.05, DIC=2200.0, ALK=2400.0)

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Show the progress of the simulation
# - Store the output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implimentation page)
# - Adapt the timestep length to reduce the run time

simulation = Simulation(model, Δt=1minutes, stop_time=100days) 

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))                       
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "data_forced"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, 
                                                        merge(model.tracers, model.auxiliary_fields), 
                                                        filename = "$filename.jld2", 
                                                        schedule = TimeInterval(1day), 
                                                        overwrite_existing = true)

simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))

simulation.callbacks[:timestep] = Callback(update_timestep!, IterationInterval(1), (c_forcing=0.2, c_adv=0.2, c_diff=0.2, w = 200/day, relaxation=0.95), TimeStepCallsite())

# ## Run!
# Finally we run the simulation
run!(simulation)

# Now we can visulise the results with some post processing to diagnose the air-sea CO₂ flux

P = FieldTimeSeries("$filename.jld2", "P")
NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
Z = FieldTimeSeries("$filename.jld2", "Z")
D = FieldTimeSeries("$filename.jld2", "D") 
DD = FieldTimeSeries("$filename.jld2", "DD")
DIC = FieldTimeSeries("$filename.jld2", "DIC")
Dᶜ = FieldTimeSeries("$filename.jld2", "Dᶜ")
DDᶜ = FieldTimeSeries("$filename.jld2", "DDᶜ")
ALK = FieldTimeSeries("$filename.jld2", "ALK")

x, y, z = nodes(P)
times = P.times

air_sea_CO₂_flux = zeros(size(P)[4])
carbon_export = zeros(size(P)[4])
for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = Boundaries.airseaflux(0.0, 0.0, t, DIC[1, 1, Nz, i], ALK[1, 1, Nz, i], t_function(1, 1, 0, t), s_function(1, 1, 0, t),  merge(Boundaries.defaults.airseaflux, (T=t_function, S=s_function, gas=:CO₂)))
    carbon_export[i] = (Dᶜ[1, 1, end-20, i]*LOBSTER.D_sinking .+ DDᶜ[1, 1, end-20, i]*LOBSTER.DD_sinking)
end

using GLMakie
f=Figure(backgroundcolor=RGBf(1, 1, 1), fontsize=30)

axP = Axis(f[1, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Phytoplankton concentration (mmol N/m³)")
hmP = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(P[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbP = Colorbar(f[1, 3], hmP)

axNO₃ = Axis(f[1, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Nitrate concentration (mmol N/m³)")
hmNO₃ = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(NO₃[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbNO₃ = Colorbar(f[1, 6], hmNO₃)

axZ = Axis(f[2, 1:2], ylabel="z (m)", xlabel="Time (days)", title="Zooplankton concentration (mmol N/m³)")
hmZ = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(Z[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbZ = Colorbar(f[2, 3], hmZ)

axD = Axis(f[2, 4:5], ylabel="z (m)", xlabel="Time (days)", title="Detritus concentration (mmol N/m³)")
hmD = GLMakie.heatmap!(times./days, float.(z[end-23:end]), float.(D[1, 1, end-23:end, 1:101])' .+ float.(DD[1, 1, end-23:end, 1:101])', interpolate=true, colormap=:batlow)
cbD = Colorbar(f[2, 6], hmD)

axfDIC = Axis(f[3, 1:4], xlabel="Time (days)", title="Air-sea CO₂ flux and Sinking", ylabel="Flux (kgCO₂/m²/year)")
hmfDIC = GLMakie.lines!(times./days, air_sea_CO₂_flux.*(12+16*2).*year/(1000*1000), label="Air-sea flux")
hmfExp = GLMakie.lines!(times./days, carbon_export.*(12+16*2).*year/(1000*1000), label="Sinking export")

f[3, 5] = Legend(f, axfDIC, "", framevisible = false)

save("examples/data_forced.png", f)

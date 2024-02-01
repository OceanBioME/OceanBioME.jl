# # One dimensional column forced by external data with carbonate chemistry
# In this example we setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution.
# This example demonstrates:
# - How to setup OceanBioME's biogeochemical models
# - How to load external forcing data
# - How to run with optional tracer sets such as carbonate chemistry
# - How to setup a non-uniform grid for better near surface resolution
# - How to visualise results

# For this example we use force by mixing layer depth and surface photosynthetically available radiation (PAR) data from the
# Mercator Ocean model and NASA VIIRS observations.

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, NetCDF, Interpolations, DataDeps, CairoMakie"
# ```

# ## Model setup
# First load the required packages
using OceanBioME
using Oceananigans, Random, Printf, NetCDF, Interpolations, DataDeps
using Oceananigans.Units
using Oceananigans.Fields: FunctionField

import Oceananigans.TurbulenceClosures: maximum_numeric_diffusivity

maximum_numeric_diffusivity(κ::NamedTuple) = maximum(maximum.(values(κ)))
maximum_numeric_diffusivity(κ::FunctionField) = maximum(κ)

const year = years = 365days # just for these idealised cases
nothing #hide

# ## Load external forcing data
# Loading the forcing data from our online copy
dd = DataDep(
    "example_data",
    "example data from subpolar re-analysis and observational products",
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
κₜ(x, y, z, t) = 2e-2 * max(1 - (z + mld_itp(mod(t, 364days)) / 2)^2 / (mld_itp(mod(t, 364days)) / 2)^2, 0) + 1e-4

nothing #hide

# ## Grid and diffusivity field
# Define the grid (in this case a non uniform grid for better resolution near the surface) and an extra Oceananigans field for the PAR to be stored in
Nz = 33
Lz = 600meters
refinement = 10 
stretching = 5.754
h(k) = (k - 1) / Nz
ζ₀(k) = 1 + (h(k) - 1) / refinement
Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))
z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

grid = RectilinearGrid(size = (1, 1, Nz), x = (0, 20meters), y = (0, 20meters), z = z_faces)

clock = Clock(; time = 0.0)

κ = FunctionField{Center, Center, Center}(κₜ, grid; clock)

# ## Biogeochemical and Oceananigans model
# Here we instantiate the LOBSTER model with carbonate chemistry and a surface flux of DIC (CO₂)

biogeochemistry = LOBSTER(; grid,
                            surface_photosynthetically_active_radiation = surface_PAR,
                            carbonates = true,
                            scale_negatives = true)

CO₂_flux = GasExchange(; gas = :CO₂)

T = FunctionField{Center, Center, Center}(t_function, grid; clock)
S = FunctionField{Center, Center, Center}(s_function, grid; clock)

model = NonhydrostaticModel(; grid, clock,
                              closure = ScalarDiffusivity(ν = κ, κ = κ),
                              biogeochemistry,
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),),
                              auxiliary_fields = (; T, S))

set!(model, P = 0.03, Z = 0.03, NO₃ = 11.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2400.0)

# ## Simulation
# Next we setup the simulation along with some callbacks that:
# - Show the progress of the simulation
# - Store the output
# - Prevent the tracers from going negative from numerical error (see discussion of this in the [positivity preservation](@ref pos-preservation) implementation page)
# - Adapt the timestep length to reduce the run time

simulation = Simulation(model, Δt = 1minutes, stop_time = 100days)

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(500))

filename = "data_forced"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, 
                                                        model.tracers,
                                                        filename = "$filename.jld2",
                                                        schedule = TimeInterval(1day),
                                                        overwrite_existing = true)

wizard = TimeStepWizard(cfl = 0.2, diffusive_cfl = 0.2,
                        max_change = 1.5, min_change = 0.75,
                        cell_advection_timescale = column_advection_timescale)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
nothing #hide

# ## Run!
# We are ready to run the simulation
run!(simulation)

# ## Load output and plot
# Now we can visualise the results with some post processing to diagnose the air-sea CO₂ flux

   P = FieldTimeSeries("$filename.jld2", "P")
 NO₃ = FieldTimeSeries("$filename.jld2", "NO₃")
   Z = FieldTimeSeries("$filename.jld2", "Z")
sPOM = FieldTimeSeries("$filename.jld2", "sPOM")
bPOM = FieldTimeSeries("$filename.jld2", "bPOM")
 DIC = FieldTimeSeries("$filename.jld2", "DIC")
 Alk = FieldTimeSeries("$filename.jld2", "Alk")

x, y, z = nodes(P)
times = P.times
nothing #hide

# We compute the  air-sea CO₂ flux at the surface (corresponding to vertical index `k = grid.Nz`) and
# the carbon export by computing how much carbon sinks below some arbitrary depth; here we use depth
# that corresponds to `k = grid.Nz - 20`.
air_sea_CO₂_flux = zeros(length(times))
carbon_export = zeros(length(times))

using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity

for (i, t) in enumerate(times)
    air_sea_CO₂_flux[i] = CO₂_flux.condition.func(0.0, 0.0, t, DIC[1, 1, grid.Nz, i], Alk[1, 1, grid.Nz, i], t_function(1, 1, 0, t), s_function(1, 1, 0, t))
    carbon_export[i] = (sPOM[1, 1, grid.Nz-20, i] * biogeochemical_drift_velocity(model.biogeochemistry, Val(:sPOM)).w[1, 1, grid.Nz-20] +
                        bPOM[1, 1, grid.Nz-20, i] * biogeochemical_drift_velocity(model.biogeochemistry, Val(:bPOM)).w[1, 1, grid.Nz-20]) * redfield(Val(:sPOM), model.biogeochemistry)
end

# Both `air_sea_CO₂_flux` and `carbon_export` are in units `mmol CO₂ / (m² s)`.

using CairoMakie

fig = Figure(size = (1000, 1500), fontsize = 20)

axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((0, times[end] / days), (-150meters, 0)))

axP = Axis(fig[1, 1]; title = "Phytoplankton concentration (mmol N/m³)", axis_kwargs...)
hmP = heatmap!(times / days, z, interior(P, 1, 1, :, :)', colormap=:batlow)
Colorbar(fig[1, 2], hmP)

axNO₃ = Axis(fig[2, 1]; title = "Nitrate concentration (mmol N/m³)", axis_kwargs...)
hmNO₃ = heatmap!(times / days, z, interior(NO₃, 1, 1, :, :)', colormap=:batlow)
Colorbar(fig[2, 2], hmNO₃)

axZ = Axis(fig[3, 1]; title = "Zooplankton concentration (mmol N/m³)", axis_kwargs...)
hmZ = heatmap!(times / days, z, interior(Z, 1, 1, :, :)', colormap=:batlow)
Colorbar(fig[3, 2], hmZ)

axD = Axis(fig[4, 1]; title = "Detritus concentration (mmol N/m³)", axis_kwargs...)
hmD = heatmap!(times / days, z, interior(sPOM, 1, 1, :, :)' .+ interior(bPOM, 1, 1, :, :)', colormap=:batlow)
Colorbar(fig[4, 2], hmD)

CO₂_molar_mass = (12 + 2 * 16) * 1e-3 # kg / mol

axfDIC = Axis(fig[5, 1], xlabel = "Time (days)", ylabel = "Flux (kgCO₂/m²/year)",
                         title = "Air-sea CO₂ flux and Sinking", limits = ((0, times[end] / days), nothing))
lines!(axfDIC, times / days, cumsum(air_sea_CO₂_flux) /1e3 * CO₂_molar_mass * year, linewidth = 3, label = "Air-sea flux")
lines!(axfDIC, times / days, cumsum(carbon_export) /1e3    * CO₂_molar_mass * year, linewidth = 3, label = "Sinking export")
Legend(fig[5, 2], axfDIC, framevisible = false)

fig

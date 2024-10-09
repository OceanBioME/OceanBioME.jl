# # [One-dimensional column example](@id OneD_column)
# In this example we setup a simple 1D column with the [LOBSTER](@ref LOBSTER) biogeochemical model and observe its evolution. The example demonstrates:
# - How to setup OceanBioME's biogeochemical models
# - How to visualise results
# This is forced by idealised mixing layer depth and surface photosynthetically available radiation (PAR) which are setup first.

# ## Install dependencies
# First we check we have the dependencies installed
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, CairoMakie"
# ```

# ## Model setup
# We load the packages and choose the default LOBSTER parameter set
using OceanBioME, Oceananigans, Printf

using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units

using OceanBioME.Sediments: sinking_flux

const year = years = 365days
nothing #hide

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic), temperaeture and euphotic layer

@inline PAR⁰(t) = 300 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 10

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))

@inline κₜ(z, t) = (1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4)

@inline temp(z, t) = 2.4 * (1 + cos(t * 2π / year + 50days)) * ifelse(z > MLD(t), 1, exp((z - MLD(t))/20)) + 8

grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (100, ), extent = (400, ))

clock = Clock(; time = 0.0)

# we can keep this in the column version where we are compleltly divorced from the physics but it be the default to compute it
# JSW will implement somewhere else and it can be pulled in and made the default at some point before merge
zₘₓₗ = FunctionField{Center, Center, Nothing}(MLD, grid; clock)
κ_field = FunctionField{Center, Center, Center}(κₜ, grid; clock)

# ## Model
# First we define the biogeochemical model including carbonate chemistry (for which we also define temperature (``T``) and salinity (``S``) fields)
# and scaling of negative tracers(see discussion in the [positivity preservation](@ref pos-preservation))
# and then setup the Oceananigans model with the boundary condition for the DIC based on the air-sea CO₂ flux.

carbon_chemistry = CarbonChemistry()

biogeochemistry = PISCES(; grid,
                           mixed_layer_depth = zₘₓₗ,
                           mean_mixed_layer_vertical_diffusivity = ConstantField(1e-2), # this is by default computed now
                           surface_photosynthetically_active_radiation = PAR⁰,
                           carbon_chemistry)

CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition(; carbon_chemistry)
O₂_flux = OxygenGasExchangeBoundaryCondition()

@info "Setting up the model..."
model = HydrostaticFreeSurfaceModel(; grid,
                                      velocities = PrescribedVelocityFields(),
                                      tracer_advection = TracerAdvection(nothing, nothing, WENOFifthOrder(grid)),
                                      momentum_advection = nothing,
                                      buoyancy = nothing,
                                      clock,
                                      closure = ScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ = κ_field),
                                      biogeochemistry,
                                      boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), O₂ = FieldBoundaryConditions(top = O₂_flux)))
                              

@info "Setting initial values..."

set!(model, P = 0.1, PChl = 0.025, PFe = 0.005,
            D = 0.01, DChl = 0.003, DFe = 0.0006, DSi = 0.004,
            Z = 0.06, M = 0.5,
            DOC = 4, 
            POC = 5.4, SFe = 0.34, 
            GOC = 8.2, BFe = 0.5, PSi = 0.04, CaCO₃ = 10^-10,
            NO₃ = 10, NH₄ = 0.1, PO₄ = 5.0, Fe = 0.6, Si = 8.6,
            DIC = 2205, Alk = 2560, O₂ = 317, S = 35)

# maybe get to 1.5hours after initial stuff
simulation = Simulation(model, Δt = 30minutes, stop_time = 5years)

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

add_callback!(simulation, progress_message, TimeInterval(10day))

# prescribe the temperature
using KernelAbstractions: @index, @kernel
using Oceananigans.Architectures: architecture
using Oceananigans.Grids: znode, Center
using Oceananigans.Utils: launch!

@kernel function fill_T!(T, grid, temp, t)
    i, j, k = @index(Global, NTuple)

    z = znode(i, j, k, grid, Center(), Center(), Center())

    @inbounds T[i, j, k] = temp(z, t)

end
function update_temperature!(simulation)
    t = time(simulation)

    launch!(architecture(grid), grid, :xyz, fill_T!, model.tracers.T, grid, temp, t)

    return nothing
end

add_callback!(simulation, update_temperature!, IterationInterval(1))

filename = "column"
simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "$filename.jld2",
                                                       schedule = TimeInterval(3day),
                                                       overwrite_existing = true)

PAR = Field(Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(biogeochemistry.light_attenuation).PAR)

internal_fields = (; biogeochemistry.underlying_biogeochemistry.calcite_saturation, 
                     biogeochemistry.underlying_biogeochemistry.euphotic_depth, 
                     PAR
                     )#biogeochemistry.underlying_biogeochemistry.mean_mixed_layer_vertical_diffusivity)

simulation.output_writers[:internals] = JLD2OutputWriter(model, internal_fields,
                                                       filename = "$(filename)_internal_fields.jld2",
                                                       schedule = TimeInterval(3day),
                                                       overwrite_existing = true)
# ## Run!
# We are ready to run the simulation
run!(simulation)


# ## Load saved output
# Now we can load the results and do some post processing to diagnose the air-sea CO₂ flux. Hopefully, this looks different to the example without kelp!

tracers = FieldDataset("$filename.jld2")
internal_fields = FieldDataset("$(filename)_internal_fields.jld2")

x, y, z = nodes(tracers["P"])
times = tracers["P"].times

# We compute the  air-sea CO₂ flux at the surface (corresponding to vertical index `k = grid.Nz`) and
# the carbon export by computing how much carbon sinks below some arbirtrary depth; here we use depth 
# that corresponds to `k = grid.Nz - 20`.
air_sea_CO₂_flux = zeros(length(times))
carbon_export = zeros(length(times))

S = ConstantField(35)

using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using CUDA

CUDA.@allowscalar for (n, t) in enumerate(times)
    clock.time = t

    k_export = floor(Int, grid.Nz + MLD(t)/minimum_zspacing(grid))
    
    air_sea_CO₂_flux[n] = CO₂_flux.condition.func(1, 1, grid, clock, (; DIC = tracers["DIC"][n], Alk = tracers["Alk"][n], T = tracers["T"][n], S))

    POC_export = -sinking_flux(1, 1, k_export, grid, model.advection.POC.z, Val(:POC), biogeochemistry, (; POC = tracers["POC"][n]))
    GOC_export = -sinking_flux(1, 1, k_export, grid, model.advection.GOC.z, Val(:GOC), biogeochemistry, (; GOC = tracers["GOC"][n]))

    carbon_export[n] = POC_export + GOC_export
end

using CairoMakie

fig = Figure(size = (4000, 2100), fontsize = 20);

start_day = 366
end_day   = length(times)

axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((times[start_day], times[end_day]) ./ days, (-200, 0)))

for (n, name) in enumerate(keys(model.tracers))
    if !(name == :S)
        i = floor(Int, (n-1)/4) + 1
        j = mod(2 * (n-1), 8) + 1
        ax = Axis(fig[i, j]; title = "$name", axis_kwargs...)
        hm = heatmap!(ax, times[start_day:end]./days, z, interior(tracers["$name"], 1, 1, :, start_day:end_day)')
        Colorbar(fig[i, j+1], hm)
        lines!(ax, times[start_day:end_day]./days, t->MLD(t*day), linewidth = 2, color = :black, linestyle = :dash)
        lines!(ax, times[start_day:end_day]./days, interior(internal_fields["euphotic_depth"], 1, 1, 1, start_day:end_day), linewidth = 2, color = :white, linestyle = :dash)
    end
end

ax = Axis(fig[7, 3]; title = "log₁₀(Calcite saturation)", axis_kwargs...)
hm = heatmap!(ax, times[start_day:end_day]./days, z, log10.(interior(internal_fields["calcite_saturation"], 1, 1, :, start_day:end_day)'))
Colorbar(fig[7, 4], hm)

ax = Axis(fig[7, 5]; title = "log₁₀(PAR)", axis_kwargs...)
hm = heatmap!(ax, times[start_day:end_day]./days, z, log10.(interior(internal_fields["PAR"], 1, 1, :, start_day:end_day)'))
Colorbar(fig[7, 6], hm)

CO₂_molar_mass = (12 + 2 * 16) * 1e-3 # kg / mol

axDIC = Axis(fig[7, 7], xlabel = "Time (days)", ylabel = "Flux (kgCO₂/m²/year)",
                         title = "Air-sea CO₂ flux and Sinking", limits = ((times[start_day], times[end_day]) ./ days, nothing))

lines!(axDIC, times[start_day:end_day] / days, air_sea_CO₂_flux[start_day:end_day] / 1e3 * CO₂_molar_mass * year, linewidth = 3, label = "Air-sea flux")
lines!(axDIC, times[start_day:end_day] / days, carbon_export[start_day:end_day] / 1e3    * CO₂_molar_mass * year, linewidth = 3, label = "Sinking below mixed layer")
Legend(fig[7, 8], axDIC, framevisible = false)

fig
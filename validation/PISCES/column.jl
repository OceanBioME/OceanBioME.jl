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

const year = years = 365days
nothing #hide

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth 
# Setting up idealised functions for PAR and diffusivity (details here can be ignored but these are typical of the North Atlantic), temperaeture and euphotic layer

@inline PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

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

T = FunctionField{Center, Center, Center}(temp, grid; clock)
S = ConstantField(35)

@info "Setting up the model..."
model = NonhydrostaticModel(; grid,
                              clock,
                              closure = ScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ = κ_field),
                              biogeochemistry,
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), O₂ = FieldBoundaryConditions(top = O₂_flux)),
                              auxiliary_fields = (; S))
                              

@info "Setting initial values..."
set!(model, P = 6.95, D = 6.95, Z = 0.695,  M = 0.695, 
            Pᶜʰˡ = 1.671,  Dᶜʰˡ = 1.671, 
            Pᶠᵉ =7e-6 * 1e9 / 1e6 * 6.95, Dᶠᵉ = 7e-6 * 1e9 / 1e6 * 6.95, 
            Dˢⁱ = 1.162767, 
            DOC = 0.0, POC = 0.0, GOC = 0.0, 
            SFe = 7e-6 * 1e9 / 1e6 *1.256, BFe =7e-6 * 1e9 / 1e6 *1.256, 
            NO₃ = 6.202, NH₄ = 0.25*6.202, 
            PO₄ = 0.8722, Fe = 1.256, Si = 7.313, 
            CaCO₃ = 0.001,
            DIC = 2139.0, Alk = 2366.0, 
            O₂ = 237.0) #Using Copernicus Data (26.665, 14.), Calcite is not correct, but this is to see it on the graphs

# ## Simulation
# Next we setup a simulation and add some callbacks that:
# - Show the progress of the simulation
# - Store the model and particles output

simulation = Simulation(model, Δt = 50minutes, stop_time = 2years)

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

add_callback!(simulation, progress_message, TimeInterval(10day))

# prescribe the temperature
function update_temperature!(simulation)
    t = time(simulation)

    T = reshape(map(z -> temp(z, t), znodes(simulation.model.grid, Center())), (1, 1, size(grid, 3)))

    set!(simulation.model.tracers.T, T)

    return nothing
end

add_callback!(simulation, update_temperature!, IterationInterval(1))

#NaN Checker function. Could be removed to improve speed, if confident of model stability
function non_zero_fields!(model) 
    @inbounds for (idx, tracer) in enumerate(model.tracers)
        for i in 1:50
            if isnan(tracer[1,1,i])
                throw("$(keys(model.tracers)[idx]) has gone NaN")
            else
                tracer[1, 1, i] = max(0, tracer[1, 1, i])
            end
        end
        
    end
    return nothing
end

simulation.callbacks[:non_zero_fields] = Callback(non_zero_fields!, callsite = UpdateStateCallsite())
filename = "column"
simulation.output_writers[:tracers] = JLD2OutputWriter(model, model.tracers,
                                                       filename = "$filename.jld2",
                                                       schedule = TimeInterval(1day),
                                                       overwrite_existing = true)

internal_fields = (; biogeochemistry.underlying_biogeochemistry.calcite_saturation, 
                     biogeochemistry.underlying_biogeochemistry.euphotic_depth, 
                     )#biogeochemistry.underlying_biogeochemistry.mean_mixed_layer_vertical_diffusivity)

simulation.output_writers[:internals] = JLD2OutputWriter(model, internal_fields,
                                                       filename = "$(filename)_internal_fields.jld2",
                                                       schedule = TimeInterval(1day),
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

using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity

for (n, t) in enumerate(times)

    clock.time = t
    compute!(T)
    
    air_sea_CO₂_flux[n] = CO₂_flux.condition.func(1, 1, grid, clock, (; DIC = tracers["DIC"][n], Alk = tracers["Alk"][n], T, S))

    POC = interior(tracers["POC"][n], 1, 1, grid.Nz-20)[1]
    wPOC = biogeochemical_drift_velocity(model.biogeochemistry, Val(:POC)).w[1, 1, grid.Nz-20]

    GOC = interior(tracers["GOC"][n], 1, 1, grid.Nz-20)[1]
    wGOC = biogeochemical_drift_velocity(model.biogeochemistry, Val(:GOC)).w[1, 1, grid.Nz-20]

    carbon_export[n] = (POC * wPOC + GOC * wGOC)
end

using CairoMakie

fig = Figure(size = (4000, 2100), fontsize = 20);

axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)", limits = ((180, times[end] / days), (-400, 0)))

for (n, name) in enumerate(keys(model.tracers))
    i = floor(Int, (n-1)/4) + 1
    j = mod(2 * (n-1), 8) + 1
    ax = Axis(fig[i, j]; title = "$name", axis_kwargs...)
    hm = heatmap!(ax, times[180:end]./days, z, interior(tracers["$name"], 1, 1, :, 180:731)')
    Colorbar(fig[i, j+1], hm)
    lines!(ax, times[180:end]./days, t->MLD(t*day), linewidth = 2, color = :black, linestyle = :dash)
    lines!(ax, times[180:end]./days, interior(internal_fields["euphotic_depth"], 1, 1, 1, 180:731), linewidth = 2, color = :white, linestyle = :dash)
end

ax = Axis(fig[7, 3]; title = "log10(Calcite saturation), looks temperature dominated", axis_kwargs...)
hm = heatmap!(ax, times[180:end]./days, z, log10.(interior(internal_fields["calcite_saturation"], 1, 1, :, 180:731)'), colorrange = (-173, -95))
Colorbar(fig[7, 4], hm)

CO₂_molar_mass = (12 + 2 * 16) * 1e-3 # kg / mol

axDIC = Axis(fig[7, 5], xlabel = "Time (days)", ylabel = "Flux (kgCO₂/m²/year)",
                         title = "Air-sea CO₂ flux and Sinking", limits = ((0, times[end] / days), nothing))

lines!(axDIC, times[180:731] / days, air_sea_CO₂_flux[180:731] / 1e3 * CO₂_molar_mass * year, linewidth = 3, label = "Air-sea flux")
lines!(axDIC, times[180:731] / days, carbon_export[180:731] / 1e3    * CO₂_molar_mass * year, linewidth = 3, label = "Sinking export")
Legend(fig[7, 6], axDIC, framevisible = false)

fig

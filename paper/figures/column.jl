# This script is not continuously tested so may be out of date
# It is compatible with the version of OceanBioME available at JOSS paper release (v0.X.Y)
using OceanBioME, Oceananigans, Printf, JLD2
using OceanBioME.SLatissimaModel: SLatissima
using Oceananigans.Units

const year = years = 365days

@inline PAR⁰(x, y, t) = 60 * (1 - cos((t + 15days) * 2π / year))*(1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days) ^ 2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)

@inline fmld1(t) = H(t, 50days, year) * (1 / (1 +exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))

@inline MLD(t) = - (10 + 340 * (1 - fmld1(year - eps(year)) * exp(- mod(t, year) / 25days) - fmld1(mod(t, year))))

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4

@inline t_function(x, y, z, t) = 2.4 * cos(t * 2π / year + 50day) + 10

# ## Grid and PAR field
# Define the grid and an extra Oceananigans field for the PAR to be stored in
Lx, Ly = 20, 20
grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200)) 

#####
##### Initial warmup
#####

CO₂_flux = GasExchange(; gas = :CO₂, temperature = t_function, salinity = (args...) -> 35)

N_fixation(x, y, t) = - 10 / year
N_upwelling = Relaxation(rate = 1/10days, mask = (x, y, z) -> ifelse(z < -150, 1.0, 0.0), target = 1.0)

model = NonhydrostaticModel(; grid,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                              biogeochemistry = LOBSTER(; grid,
                                                          surface_phytosynthetically_active_radiation = PAR⁰,
                                                          carbonates = true,
                                                          variable_redfield = true,
                                                          advection_schemes = (sPOM = WENO(grid), bPOM = WENO(grid))),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),
                                                     NO₃ = FieldBoundaryConditions(top = FluxBoundaryCondition(N_fixation))),
                              forcing = (NO₃ = N_upwelling, ),
                              advection = nothing)

set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

simulation = Simulation(model, Δt = 10minutes, stop_time = 3year) 

simulation.callbacks[:nan_checker] = Callback(Oceananigans.Simulations.NaNChecker(; fields = model.tracers, erroring = true))

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time)) 

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())

wizard = TimeStepWizard(cfl = 0.15, diffusive_cfl = 0.15, max_change = 2.0, min_change = 0.5, cell_diffusion_timescale = column_diffusion_timescale, cell_advection_timescale = column_advection_timescale)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=SpecifiedTimes([3years, 5years - 30days]),
                                                        prefix = "warmup_checkpoint", overwrite_existing = true)

run!(simulation) # initial warmup

simulation.stop_time = 5years - 30days

filename = "column"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)

run!(simulation, pickup = false)

# ## Kelp Particle setup
@info "Setting up kelp particles"
n = 5 # number of kelp fronds
z₀ = [-21:5:-1;] * 1.0 # depth of kelp fronds

particles = SLatissima(; x = ones(n) * Lx / 2, y = ones(n) * Ly / 2, z = z₀, 
                         A = ones(n) * 10.0, N = ones(n) * 0.014, C =  ones(n) * 0.4, 
                         latitude = 57.5,
                         scalefactor = 50.0, 
                         pescribed_temperature = t_function)

CO₂_flux = GasExchange(; gas = :CO₂, temperature = t_function, salinity = (args...) -> 35)

N_fixation(x, y, t) = - 10 / year
N_upwelling = Relaxation(rate = 1/10days, mask = (x, y, z) -> ifelse(z < -150, 1.0, 0.0), target = 1.0)

model = NonhydrostaticModel(; grid,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                              biogeochemistry = LOBSTER(; grid,
                                                          surface_phytosynthetically_active_radiation = PAR⁰,
                                                          carbonates = true,
                                                          variable_redfield = true,
                                                          particles,
                                                          advection_schemes = (sPOM = WENO(grid), bPOM = WENO(grid))),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),
                                                     NO₃ = FieldBoundaryConditions(top = FluxBoundaryCondition(N_fixation))),
                              forcing = (NO₃ = N_upwelling, ),
                              advection = nothing)

model.clock.time = 5year - 30days

set!(model, "warmup_checkpoint_iteration634986.jld2")

simulation = Simulation(model, Δt=1minutes, stop_time=7years) 

simulation.callbacks[:nan_checker] = Callback(Oceananigans.Simulations.NaNChecker(; fields = model.tracers, erroring = true))

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                                        iteration(sim),
                                                        prettytime(sim),
                                                        prettytime(sim.Δt),
                                                        prettytime(sim.run_wall_time))  

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "kelp"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, merge(model.tracers, model.auxiliary_fields), filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)
simulation.output_writers[:particles] = JLD2OutputWriter(model, (; particles), filename = "$(filename)_particles.jld2", schedule = TimeInterval(1day), overwrite_existing=true)

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# ## Run!
# Finally we run the simulation
run!(simulation)

using OceanBioME, Oceananigans, Printf
using Oceananigans.Units

year = years = 365days # just for these idealised cases

@inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - 70)/10)) / 2 + 1e-4

Lx, Ly = 20, 20
grid = RectilinearGrid(size=(1, 1, 50), extent=(Lx, Ly, 200)) 

# Specify the boundary conditions for DIC and O₂ based on the air-sea CO₂ and O₂ flux
CO₂_flux = GasExchange(; gas = :CO₂, temperature = (args...) -> 12, salinity = (args...) -> 35)

# this rate doesn't seem excessive compared to the nitrogen fixation rates from https://www.nature.com/articles/nature05392
# riverine input is presumably more significant than this in coastal waters
N_fixation(x, y, t) = - 10 / year
N_upwelling = Relaxation(rate = 1/10days, mask = (x, y, z) -> ifelse(z < -150, 1.0, 0.0), target = 1.0)

# alkalinity - restoring to a north atlantic climatology otherwise we depleat it
Alk_forcing = Relaxation(rate = 1/day, target = 2409.0)

model = NonhydrostaticModel(; grid,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                              biogeochemistry = LOBSTER(; grid,
                                                          surface_phytosynthetically_active_radiation = (x, y, t) -> 100,
                                                          carbonates = true,
                                                          open_bottom = true,),
                              boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),
                                                     NO₃ = FieldBoundaryConditions(top = FluxBoundaryCondition(N_fixation))),
                              forcing = (NO₃ = N_upwelling, Alk = Alk_forcing),
                              advection = nothing)

set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2400.0)

simulation = Simulation(model, Δt=5minutes, stop_time=10year) 

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))      
                                                                  
simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(100))

filename = "steady_state"
simulation.output_writers[:profiles] = JLD2OutputWriter(model, model.tracers, filename = "$filename.jld2", schedule = TimeInterval(1day), overwrite_existing=true)

scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())

run!(simulation)

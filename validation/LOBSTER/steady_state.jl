using Oceananigans, OceanBioME, Printf
using Oceananigans.Units

grid = RectilinearGrid(; topology = (Flat, Flat, Bounded), size = (1, ), extent = (40, ))

# biogeochemistry
biogeochemistry = LOBSTER(; grid, 
                            carbonates = true, oxygen = true, variable_redfield = true, 
                            open_bottom = true,
                            surface_phytosynthetically_active_radiation = (x, y, t) -> 80)

DIC_bcs = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂))
O₂_bcs = FieldBoundaryConditions(top = GasExchange(; gas = :O₂))

# nitrate - restoring to a (made up) climatology otherwise we depleat Nutrients
NO₃_forcing = Relaxation(rate = 1/10days, target = 3.0)

# alkalinity - restoring to a (made up) climatology otherwise we depleat it
Alk_forcing = Relaxation(rate = 1/day, target = 2409.0)

# define model
model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              forcing = (NO₃ = NO₃_forcing, Alk = Alk_forcing),
                              boundary_conditions = (DIC = DIC_bcs, O₂ = O₂_bcs),
                              timestepper = :RungeKutta3,
                              tracers = (:T, :S))

set!(model, P = 0.03, Z = 0.03, NO₃ = 11.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2400.0, O₂ = 240.0, T=20, S=35)

simulation = Simulation(model, Δt=1000.0, stop_time=10years)

# Print a progress message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, max(|w|) = %.1e ms⁻¹, wall time: %s\n",
                                iteration(sim), prettytime(sim), prettytime(sim.Δt),
                                maximum(abs, sim.model.velocities.w), prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, IterationInterval(50))

filename = "steady_state"

simulation.output_writers[:tracers] =
    JLD2OutputWriter(model, model.tracers,
                     filename = filename * ".jld2",
                     schedule = TimeInterval(0.5days),
                     overwrite_existing = true)

run!(simulation)

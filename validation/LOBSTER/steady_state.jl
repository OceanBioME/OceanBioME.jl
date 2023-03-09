using Oceananigans, OceanBioME, Printf
using Oceananigans.Units

grid = RectilinearGrid(; topology = (Flat, Flat, Bounded), size = (16, ), extent = (100, ))

# biogeochemistry
biogeochemistry = LOBSTER(; grid, 
                            carbonates = true,
                            open_bottom = true,
                            surface_phytosynthetically_active_radiation = (x, y, t) -> 100)

DIC_bcs = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂))

# nitrate - restoring to a (made up) climatology otherwise we depleat Nutrients
NO₃_forcing = Relaxation(rate = 1/10days, target = 4.0)

# alkalinity - restoring to a north atlantic climatology otherwise we depleat it
Alk_forcing = Relaxation(rate = 1/day, target = 2409.0)

# define model
κₜ(x, y, z, t) = 1e-4 + 1e-2 * (1 + tanh((z + 70) / 70)) / (1 + tanh(1))

model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              forcing = (NO₃ = NO₃_forcing, Alk = Alk_forcing),
                              boundary_conditions = (DIC = DIC_bcs, ),
                              timestepper = :RungeKutta3,
                              tracers = (:T, :S),
                              advection = nothing,
                              closure = ScalarDiffusivity(ν = κₜ, κ = κₜ))

set!(model, P = 0.4686, Z = 0.5363, 
            NO₃ = 2.3103, NH₄ = 0.0010, 
            DIC = 2106.9, Alk = 2408.9, 
            DOM = 0.8115,
            sPOM = 0.2299,
            bPOM = 0.0103)

simulation = Simulation(model, Δt=500.0, stop_time=2years)

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

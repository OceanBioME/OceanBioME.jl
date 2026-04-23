using OceanBioME, Printf, Oceananigans, Oceananigans.Units

using Oceananigans.Architectures: on_architecture

using Oceananigans.Solvers: FFTBasedPoissonSolver

using Random

Random.seed!(43)

h(k) = (k - 1) / Nz

ζ₀(k) = 1 + (h(k) - 1) / refinement

Σ(k) = (1 - exp(-stretching * h(k))) / (1 - exp(-stretching))

z_faces(k) = Lz * (ζ₀(k) * Σ(k) - 1)

Nx, Ny, Nz = 512, 512, 64
Lx, Ly, Lz = 1kilometer, 1kilometer, 140.0

refinement = 1.8
stretching = 3

arch = GPU()

grid = RectilinearGrid(arch; size = (Nx, Ny, Nz), x = (0, Lx), y = (0, Ly), z = z_faces)

coriolis = FPlane(f = 1e-4) # [s⁻¹]

background_state_parameters = ( M = 1e-4,       # s⁻¹, geostrophic shear
                                f = coriolis.f, # s⁻¹, Coriolis parameter
                                N = 1e-4,       # s⁻¹, buoyancy frequency
                                H = grid.Lz )

# We assume a background buoyancy ``B`` with a constant stratification and also a constant lateral
# gradient (in the zonal direction). The background velocity components ``U`` and ``V`` are prescribed
# so that the thermal wind relationship is satisfied, that is, ``f \partial_z U = - \partial_y B`` and
# ``f \partial_z V = \partial_x B``.
B(x, y, z, t, p) = p.M ^ 2 * x + p.N ^ 2 * (z + p.H)
V(x, y, z, t, p) = p.M ^ 2 / p.f * (z + p.H)

V_field = BackgroundField(V, parameters = background_state_parameters)
B_field = BackgroundField(B, parameters = background_state_parameters)

νᵥ = κᵥ = 1e-4 # [m² s⁻¹]
closure = ScalarDiffusivity(ν = νᵥ, κ = κᵥ)

@info "Setting up kelp particles"

n = 36 # must be a square number

x = on_architecture(arch, [repeat([Lx / (sqrt(n) + 1) * n for n in 1:Int(sqrt(n))], 1, Int(sqrt(n)))...])
y = on_architecture(arch, [repeat([Ly / (sqrt(n) + 1) * n for n in 1:Int(sqrt(n))], 1, Int(sqrt(n)))'...])
z = on_architecture(arch, zeros(Float64, n))

particles = SLatissima(; architecture = arch,
                         x, y, z,
                         A = on_architecture(arch, 5.0 .* ones(n)), N = on_architecture(arch, 0.01 .* ones(n)), C = on_architecture(arch, 0.18 .* ones(n)),
                         latitude = 43.3,
                         scalefactor = 500.0,
                         pescribed_temperature = (args...) -> 12.0)

biogeochemistry = LOBSTER(; grid,
                            carbonates = true,
                            open_bottom = true,
                            particles,
                            scale_negatives = true,
			                surface_photosynthetically_active_radiation = (x, y, t) -> 100)

DIC_bcs = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂, temperature = (args...) -> 12, salinity = (args...) -> 35))

# Model instantiation
model = NonhydrostaticModel(grid;
                            biogeochemistry,
                            boundary_conditions = (DIC = DIC_bcs, ),
                            advection = WENO(grid),
                            timestepper = :RungeKutta3,
                            coriolis,
                            tracers = :b,
                            buoyancy = BuoyancyTracer(),
                            background_fields = (b = B_field, v = V_field),
                            closure)

model.clock.time = 50days

Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v=vᵢ, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2409.0)

Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())

Δt₀ = 0.6 * min(Δx, Δy, Δz) / V(0, 0, 0, 0, background_state_parameters)

simulation = Simulation(model, Δt = Δt₀, stop_time = 60days)

# Adapt the time step while keeping the CFL number fixed.
wizard = TimeStepWizard(cfl = 0.6, diffusive_cfl = 0.6, max_Δt = 30minutes, min_change = 0.1, max_change = 2)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))

# Create a progress message.
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(sim.run_wall_time),
                        prettytime(sim.Δt))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

u, v, w = model.velocities

# Periodically save the velocities and vorticity to a file.
simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.tracers, (; u, v, w));
                                                      schedule = TimeInterval(6hours),
                                                      filename = "eady.jld2",
                                                      overwrite_existing = true)

simulation.output_writers[:particles] = JLD2OutputWriter(model, (; particles);
                                                         schedule = TimeInterval(6hours),
                                                         filename = "eady_particles.jld2",
                                                         overwrite_existing = true)

run!(simulation)

simulation.output_writers[:checkpointer] = Checkpointer(model, schedule=IterationInterval(1), prefix="model_checkpoint")

run!(simulation)

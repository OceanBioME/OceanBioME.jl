include("dependencies_for_runtests.jl")

using Oceananigans.Architectures: on_architecture

# total tracer in a grid made of 1 m³ cells (mmol)
total_tracer(tracer) = sum(on_architecture(CPU(), interior(tracer)))

particle_radii(particles) = on_architecture(CPU(), particles.fields.r)

function dissolving_particle_model(grid, particles)
    biogeochemistry = NutrientsPlanktonDetritus(grid; inorganic_carbon = CarbonateSystem(), particles)

    model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing)

    set!(model, DIC = 2000, Alk = 2300)

    return model
end

@testset "Dissolving particle dissolution" begin
    grid = RectilinearGrid(architecture; size = (1, 1, 1), extent = (1, 1, 1))

    k  = 1e-6 # m/s
    r₀ = 1e-4 # m
    Δt = 0.1  # so the particle dissolves in 1000 timesteps
    scalefactor = 1e6

    particles = DissolvingParticles(1; grid,
                                    advection = nothing,
                                    scalefactors = [scalefactor],
                                    mineral_parameters = (; dissolution_rate = k,
                                                            released_tracers = (; Alk = 2, DIC = 1)))

    set!(particles, r = r₀, x = 0.5, y = 0.5, z = -0.5)

    model = dissolving_particle_model(grid, particles)

    mineral = particles.biogeochemistry
    initial_mineral = scalefactor * 4/3 * π * r₀^3 * mineral.density / mineral.molar_mass * 1000 # mmol

    # the radius shrinks at the dissolution rate
    for _ in 1:500
        time_step!(model, Δt)
    end

    @test particle_radii(particles)[1] ≈ r₀ - 500 * k * Δt

    # and stops at zero rather than overshooting
    for _ in 1:600
        time_step!(model, Δt)
    end

    @test particle_radii(particles)[1] == 0

    # all of the mineral is released, in the right ratio, up to the timestepping error
    released_Alk = total_tracer(model.tracers.Alk) - 2300
    released_DIC = total_tracer(model.tracers.DIC) - 2000

    @test isapprox(released_Alk, 2 * initial_mineral; rtol = 10 * k * Δt / r₀)
    @test released_Alk ≈ 2 * released_DIC

    # and nothing more is released once it has dissolved
    time_step!(model, Δt)

    @test total_tracer(model.tracers.Alk) - 2300 == released_Alk
end

@testset "Dissolving particle sinking" begin
    grid = RectilinearGrid(architecture; size = (1, 1, 10), extent = (2, 2, 10))

    w = 0.01 # m/s

    particles = DissolvingParticles(2; grid, mineral_parameters = (; sinking_speed = w, dissolution_rate = 0))

    # the second particle has already dissolved so shouldn't sink
    set!(particles, r = [1e-4, 0], x = 1, y = 1, z = -1)

    model = dissolving_particle_model(grid, particles)

    for _ in 1:10
        time_step!(model, 10)
    end

    z = on_architecture(CPU(), particles.z)

    @test z[1] ≈ -1 - 10 * 10 * w
    @test z[2] ≈ -1

    stokes = StokesSinking(; particle_density = 2025, water_density = 1025,
                             dynamic_viscosity = 1e-3, gravitational_acceleration = 9)

    @test stokes(1e-4) ≈ 2 * 1000 * 9 * 1e-8 / (9 * 1e-3)
end

include("dependencies_for_runtests.jl")

using Oceananigans.Architectures: on_architecture
using OceanBioME.Particles: BiogeochemicalParticles, required_particle_fields, coupled_tracers

# total tracer in a grid made of 1 m³ cells (mmol)
total_tracer(tracer) = sum(on_architecture(CPU(), interior(tracer)))

# total of a released tracer held in the particles (mmol)
function particle_content(particles, tracer)
    mineral = particles.biogeochemistry
    r = on_architecture(CPU(), particles.fields.r)
    scalefactors = on_architecture(CPU(), particles.scalefactors)

    molar_volume = mineral.molar_mass / mineral.density
    ν = mineral.released_tracers[tracer]

    return sum(@. scalefactors * 4/3 * π * r^3 / molar_volume * ν * 1000)
end

function dissolving_particle_model(grid, particles)
    biogeochemistry = NutrientsPlanktonDetritus(grid; inorganic_carbon = CarbonateSystem(), particles)

    model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing)

    set!(model, DIC = 2000, Alk = 2300)

    return model
end

@testset "Dissolving particle setup" begin
    grid = RectilinearGrid(architecture; size = (2, 2, 2), extent = (2, 2, 2))

    particles = DissolvingParticles(3; grid)

    @test particles isa BiogeochemicalParticles
    @test particles.biogeochemistry isa DissolvingMineral
    @test length(particles) == 3
    @test required_particle_fields(particles) == (:r, )
    @test coupled_tracers(particles) == (:Alk, )

    mineral = DissolvingMineral(Float32; released_tracers = (; Alk = 2, DIC = 1), dissolution_rate = r -> 1e-12 / r)

    @test mineral.density isa Float32
    @test mineral.released_tracers === (Alk = 2f0, DIC = 1f0)
    @test mineral.sinking_speed isa StokesSinking{Float32}
    @test coupled_tracers(mineral) == (:Alk, :DIC)

    stokes = StokesSinking(; particle_density = 2025, water_density = 1025,
                             dynamic_viscosity = 1e-3, gravitational_acceleration = 9)

    @test stokes(1e-4) ≈ 2 * 1000 * 9 * 1e-8 / (9 * 1e-3)
end

@testset "Dissolving particle dissolution and conservation" begin
    grid = RectilinearGrid(architecture; size = (2, 2, 2), extent = (2, 2, 2))

    k  = 1e-8 # m/s
    r₀ = 1e-4 # m
    Δt = 1

    particles = DissolvingParticles(2; grid,
                                    advection = nothing,
                                    scalefactors = [1e6, 2e6],
                                    mineral_parameters = (; dissolution_rate = k,
                                                            released_tracers = (; Alk = 2, DIC = 1)))

    set!(particles, r = r₀, x = [0.5, 1.5], y = 0.5, z = -0.5)

    model = dissolving_particle_model(grid, particles)

    initial_Alk = total_tracer(model.tracers.Alk) + particle_content(particles, :Alk)
    initial_DIC = total_tracer(model.tracers.DIC) + particle_content(particles, :DIC)
    initial_DIC_in_water = total_tracer(model.tracers.DIC)
    initial_Alk_in_water = total_tracer(model.tracers.Alk)

    for _ in 1:100
        time_step!(model, Δt)
    end

    # the radius shrinks at the (constant) dissolution rate
    @test all(on_architecture(CPU(), particles.fields.r) .≈ r₀ - 100 * k * Δt)

    # released into the cells containing the particles in proportion to their scalefactors
    Alk = on_architecture(CPU(), interior(model.tracers.Alk))
    @test Alk[1, 1, 2] > 2300
    @test (Alk[2, 1, 2] - 2300) ≈ 2 * (Alk[1, 1, 2] - 2300)
    @test Alk[1, 1, 1] == 2300

    # released in the right ratio
    released_Alk = total_tracer(model.tracers.Alk) - initial_Alk_in_water
    released_DIC = total_tracer(model.tracers.DIC) - initial_DIC_in_water
    @test released_Alk ≈ 2 * released_DIC

    # conserved up to the timestepping error, which is of order kΔt / r
    final_Alk = total_tracer(model.tracers.Alk) + particle_content(particles, :Alk)
    final_DIC = total_tracer(model.tracers.DIC) + particle_content(particles, :DIC)

    @test isapprox(final_Alk - initial_Alk, 0; atol = 1e-3 * released_Alk)
    @test isapprox(final_DIC - initial_DIC, 0; atol = 1e-3 * released_DIC)
end

@testset "Dissolving particles dissolve completely" begin
    grid = RectilinearGrid(architecture; size = (2, 2, 2), extent = (2, 2, 2))

    k  = 1e-6 # m/s
    r₀ = 1e-4 # m
    Δt = 0.1  # so the particles dissolve in 1000 timesteps

    particles = DissolvingParticles(1; grid,
                                    advection = nothing,
                                    scalefactors = [1e6],
                                    mineral_parameters = (; dissolution_rate = k))

    set!(particles, r = r₀, x = 0.5, y = 0.5, z = -0.5)

    model = dissolving_particle_model(grid, particles)

    initial_content = particle_content(particles, :Alk)

    for _ in 1:1100
        time_step!(model, Δt)
    end

    # the radius stops at zero rather than overshooting
    @test all(on_architecture(CPU(), particles.fields.r) .== 0)

    released_Alk = total_tracer(model.tracers.Alk) - 2300 * 8

    @test isapprox(released_Alk, initial_content; rtol = 10 * k * Δt / r₀)

    # and nothing more is released once they've dissolved
    time_step!(model, Δt)

    @test total_tracer(model.tracers.Alk) - 2300 * 8 == released_Alk
end

@testset "Dissolving particle sinking" begin
    grid = RectilinearGrid(architecture; size = (2, 2, 10), extent = (2, 2, 10))

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

    # particles stay in the domain when they reach the bottom
    for _ in 1:100
        time_step!(model, 10)
    end

    z = on_architecture(CPU(), particles.z)

    @test -10 <= z[1] < -9.9

    # by default they sink by Stokes' law
    particles = DissolvingParticles(1; grid, mineral_parameters = (; dissolution_rate = 0))

    set!(particles, r = 5e-5, x = 1, y = 1, z = -1)

    model = dissolving_particle_model(grid, particles)

    time_step!(model, 10)

    z = on_architecture(CPU(), particles.z)

    @test z[1] ≈ -1 - 10 * particles.biogeochemistry.sinking_speed(5e-5)
end

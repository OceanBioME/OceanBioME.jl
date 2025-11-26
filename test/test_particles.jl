include("dependencies_for_runtests.jl")

using OceanBioME.Particles: BiogeochemicalParticles

import OceanBioME.Particles: required_particle_fields, 
                             required_tracers, 
                             coupled_tracers

import OceanBioME: required_biogeochemical_tracers,
                   required_biogeochemical_auxiliary_fields

struct SimpleParticleBiogeochemistry end

@inline (::SimpleParticleBiogeochemistry)(::Val{:A}, t, A) = one(t)

@inline required_particle_fields(::SimpleParticleBiogeochemistry) = (:A, )

@inline required_tracers(::SimpleParticleBiogeochemistry) = tuple()
@inline coupled_tracers(::SimpleParticleBiogeochemistry) = tuple()

struct NothingBiogeochemistry end
(::NothingBiogeochemistry)(args...) = 0
required_biogeochemical_tracers(::NothingBiogeochemistry) = ()
required_biogeochemical_auxiliary_fields(::NothingBiogeochemistry) = ()

grid = RectilinearGrid(architecture; size = (3, 3, 3), extent = (3, 3, 3))

@testset "Testing basic particle property integration" begin
    #for TX in (Bounded, Flat) # to check it works with different topologies

    particle_biogeochemistry = SimpleParticleBiogeochemistry()

    particles = BiogeochemicalParticles(3; grid, biogeochemistry = particle_biogeochemistry, advection = nothing)

    @test length(particles) == 3

    biogeochemistry = Biogeochemistry(NothingBiogeochemistry(); particles)

    model = NonhydrostaticModel(; grid, biogeochemistry)
    
    time_step!(model, 1)

    @test all(particles.fields.A .== 1)

    set!(particles, x = 1, A = fill(5, 3))

    @test all(particles.fields.A .≈ 5) & all(particles.x .== 1)
end

@inline required_tracers(::SimpleParticleBiogeochemistry) = (:B, )

@inline (::SimpleParticleBiogeochemistry)(::Val{:A}, t, A, B) = 0.1 * B

@testset "Testing particle tracer detection" begin
    particle_biogeochemistry = SimpleParticleBiogeochemistry()

    particles = BiogeochemicalParticles(3; grid, biogeochemistry = particle_biogeochemistry, advection = nothing)

    biogeochemistry = Biogeochemistry(NothingBiogeochemistry(); particles)

    model = NonhydrostaticModel(; grid, biogeochemistry, tracers = :B)

    set!(model, B = 1)
    
    time_step!(model, 1)

    @test all(particles.fields.A .≈ 0.1)

    # also these shouldn't have moved
    @test all(particles.x .== 0)
end

@testset "Testing particle lagrangian advection" begin
    particle_biogeochemistry = SimpleParticleBiogeochemistry()

    particles = BiogeochemicalParticles(3; grid, biogeochemistry = particle_biogeochemistry)

    biogeochemistry = Biogeochemistry(NothingBiogeochemistry(); particles)

    model = NonhydrostaticModel(; grid, biogeochemistry, tracers = :B)

    set!(model, B = 1, u = 0.1, v = 0.2)
    
    time_step!(model, 1)

    @test all(particles.fields.A .≈ 0.1)
    @test all(particles.x .≈ 0.1) && all(particles.y .≈ 0.2) && all(particles.z .≈ 0)
end

coupled_tracers(::SimpleParticleBiogeochemistry) = (:B, )

@inline (::SimpleParticleBiogeochemistry)(::Val{:A}, t, A, B) = 0.1
@inline (::SimpleParticleBiogeochemistry)(::Val{:B}, t, A, B) = -0.1

@testset "Testing particle-tracer uptake" begin
    particle_biogeochemistry = SimpleParticleBiogeochemistry()

    particles = BiogeochemicalParticles(3; grid, biogeochemistry = particle_biogeochemistry, advection = nothing)

    biogeochemistry = Biogeochemistry(NothingBiogeochemistry(); particles)

    model = NonhydrostaticModel(; grid, biogeochemistry, tracers = :B)

    set!(model, B = 1)

    # since the 0, 0, 0 point is ambigously closest to both 3, 3, 3 and 1, 1, 3 (and the logic makes it go to 3, 3, 3)
    set!(particles, x = 0.5, y = 0.5)

    time_step!(model, 1)

    @test all(particles.fields.A .≈ 0.1)

    @test interior(model.tracers.B, 1, 1, 3) .≈ 0.7 # 1 - 3 * 0.1
end


@testset "Testing particle-tracer uptake with scalefactors" begin
    particle_biogeochemistry = SimpleParticleBiogeochemistry()

    particles = BiogeochemicalParticles(3; grid, 
                                           biogeochemistry = particle_biogeochemistry, 
                                           advection = nothing,
                                           scalefactors = [1, 0.5, 0.5])

    biogeochemistry = Biogeochemistry(NothingBiogeochemistry(); particles)

    model = NonhydrostaticModel(; grid, biogeochemistry, tracers = :B)

    set!(model, B = 1)

    # since the 0, 0, 0 point is ambigously closest to both 3, 3, 3 and 1, 1, 3 (and the logic makes it go to 3, 3, 3)
    set!(particles, x = 0.5, y = 0.5)
    
    time_step!(model, 1)

    @test all(particles.fields.A .≈ 0.1)

    @test interior(model.tracers.B, 1, 1, 3) .≈ 0.8
end
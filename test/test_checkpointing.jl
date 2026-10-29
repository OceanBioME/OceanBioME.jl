include("dependencies_for_runtests.jl")

using Oceananigans: prognostic_state, prognostic_fields
using Oceananigans.Fields: AbstractField, FunctionField
using Oceananigans.Units

using OceanBioME.Particles: BiogeochemicalParticles

import OceanBioME.Particles: required_particle_fields,
                             required_tracers,
                             coupled_tracers

import OceanBioME: required_biogeochemical_tracers,
                   required_biogeochemical_auxiliary_fields

"""
    direct_and_picked_up_simulations(build_simulation, N)

Return a simulation that is `run!` straight to iteration `2N` and one that is picked up from
a checkpoint written at iteration `N` and then `run!` to iteration `2N`.
"""
function direct_and_picked_up_simulations(build_simulation, N)
    direct = build_simulation(2N)
    run!(direct)

    dir = mktempdir()

    first_half = build_simulation(N)
    first_half.output_writers[:checkpointer] = Checkpointer(first_half.model; schedule = IterationInterval(N), dir)
    run!(first_half)

    picked_up = build_simulation(2N)
    picked_up.output_writers[:checkpointer] = Checkpointer(picked_up.model; schedule = IterationInterval(N), dir)
    set!(picked_up; checkpoint = :latest)
    run!(picked_up)

    rm(dir; recursive = true)

    return direct, picked_up
end

on_host(field::AbstractField) = Array(interior(field))
on_host(array) = Array(array)

function test_matching(a::NamedTuple, b::NamedTuple)
    @test keys(a) == keys(b)

    for name in keys(a)
        @test all(isapprox.(on_host(a[name]), on_host(b[name])))
    end

    return nothing
end

function test_matching_models(a, b)
    @test a.clock.iteration == b.clock.iteration
    @test a.clock.time == b.clock.time

    test_matching(prognostic_fields(a), prognostic_fields(b))
    test_matching(a.timestepper.Gⁿ, b.timestepper.Gⁿ)
    hasproperty(a.timestepper, :G⁻) && test_matching(a.timestepper.G⁻, b.timestepper.G⁻)

    return nothing
end

#####
##### Sediments
#####

sediment_initial_conditions = (N = 10, NO₃ = 10, NH₄ = 1, P = 0.1, Z = 0.05, D = 1, sPOM = 1, bPOM = 1, O₂ = 300, T = 15)

function sediment_simulation(stop_iteration; model_type, timestepper, sediment_biogeochemistry)
    grid = RectilinearGrid(architecture; size = (3, 3, 16), extent = (10, 10, 160))

    biogeochemistry = sediment_biogeochemistry(grid)

    model = model_type(grid; biogeochemistry, timestepper, buoyancy = nothing, tracers = ())

    initial_conditions = (name => value for (name, value) in pairs(sediment_initial_conditions) if name in keys(model.tracers))

    set!(model; initial_conditions...)

    return Simulation(model; Δt = 20minutes, stop_iteration, verbose = false)
end

npzd_instant_remineralisation(grid) = NPZD(grid; sediment = InstantRemineralisationSediment(grid))
lobster_multi_g(grid) = LOBSTER(grid; sediment = SimpleMultiGSediment(grid), oxygen = Oxygen())

# The split Runge-Kutta time stepper isn't supported by sediments, so they are stepped independently
sediment_model_configurations = ((NonhydrostaticModel, :RungeKutta3),
                                 (HydrostaticFreeSurfaceModel, :QuasiAdamsBashforth2),
                                 (HydrostaticFreeSurfaceModel, :SplitRungeKutta3))

@testset "Checkpointing sediments" begin
    for (model_type, timestepper) in sediment_model_configurations,
        sediment_biogeochemistry in (npzd_instant_remineralisation, lobster_multi_g)

        @info "Testing sediment checkpointing with $model_type ($timestepper) and $sediment_biogeochemistry on $architecture"

        @testset "$model_type ($timestepper), $sediment_biogeochemistry" begin
            build_simulation(stop_iteration) = sediment_simulation(stop_iteration; model_type, timestepper, sediment_biogeochemistry)

            direct, picked_up = direct_and_picked_up_simulations(build_simulation, 5)

            test_matching_models(direct.model, picked_up.model)

            direct_sediment = direct.model.biogeochemistry.sediment
            picked_up_sediment = picked_up.model.biogeochemistry.sediment

            test_matching(direct_sediment.fields, picked_up_sediment.fields)
            test_matching(direct_sediment.timestepper.Gⁿ, picked_up_sediment.timestepper.Gⁿ)
            test_matching(direct_sediment.timestepper.G⁻, picked_up_sediment.timestepper.G⁻)

            @test direct_sediment.clock.time == picked_up_sediment.clock.time == direct.model.clock.time
        end
    end
end

#####
##### Particles
#####

struct GrowingParticleBiogeochemistry end

@inline (::GrowingParticleBiogeochemistry)(::Val{:A}, t, A, B) = 0.1 * B
@inline (::GrowingParticleBiogeochemistry)(::Val{:B}, t, A, B) = - 0.1 * B

@inline required_particle_fields(::GrowingParticleBiogeochemistry) = (:A, )
@inline required_tracers(::GrowingParticleBiogeochemistry) = (:B, )
@inline coupled_tracers(::GrowingParticleBiogeochemistry) = (:B, )

# not an `AbstractBiogeochemistry`, so has no state of its own to checkpoint
struct TracerOnlyBiogeochemistry end

(::TracerOnlyBiogeochemistry)(args...) = 0
required_biogeochemical_tracers(::TracerOnlyBiogeochemistry) = ()
required_biogeochemical_auxiliary_fields(::TracerOnlyBiogeochemistry) = ()

function particle_simulation(stop_iteration)
    grid = RectilinearGrid(architecture; size = (3, 3, 3), extent = (3, 3, 3))

    particles = BiogeochemicalParticles(3; grid, biogeochemistry = GrowingParticleBiogeochemistry())

    set!(particles, x = [0.5, 1.5, 2.5], y = [2.5, 0.5, 1.5], z = [-0.5, -1.5, -2.5])

    biogeochemistry = Biogeochemistry(TracerOnlyBiogeochemistry(); particles)

    model = NonhydrostaticModel(grid; biogeochemistry, tracers = :B)

    set!(model, B = 1, u = 0.1, v = 0.2)

    return Simulation(model; Δt = 1, stop_iteration, verbose = false)
end

@testset "Checkpointing particles" begin
    @info "Testing particle checkpointing on $architecture"

    direct, picked_up = direct_and_picked_up_simulations(particle_simulation, 5)

    test_matching_models(direct.model, picked_up.model)

    direct_particles = direct.model.biogeochemistry.particles
    picked_up_particles = picked_up.model.biogeochemistry.particles

    test_matching((; direct_particles.x, direct_particles.y, direct_particles.z),
                  (; picked_up_particles.x, picked_up_particles.y, picked_up_particles.z))

    test_matching(direct_particles.fields, picked_up_particles.fields)
end

#####
##### PISCES
#####

function pisces_simulation(dir)
    grid = RectilinearGrid(architecture, topology = (Flat, Flat, Bounded), size = (4, ), extent = (100, ))

    biogeochemistry = PISCES(; grid)

    # `NonhydrostaticModel` merges the biogeochemical auxiliary fields into its own, which Oceananigans
    # can't restore when they include an `AbstractOperation` (like the `PAR` of the default light attenuation)
    model = HydrostaticFreeSurfaceModel(grid; biogeochemistry)

    simulation = Simulation(model; Δt = 1, stop_iteration = 2, verbose = false)

    simulation.output_writers[:checkpointer] = Checkpointer(model; schedule = IterationInterval(2), dir)

    return simulation
end

@testset "Checkpointing PISCES" begin
    @info "Testing PISCES checkpointing on $architecture"

    dir = mktempdir()

    simulation = pisces_simulation(dir)

    pisces = simulation.model.biogeochemistry.underlying_biogeochemistry

    # the mixed layer depth isn't computed by PISCES, so it may have been set by the user (e.g. in a callback)
    set!(pisces.mixed_layer_depth, -25)

    # the silicate climatology defaults to a `ConstantField` which has nothing to checkpoint
    @test isnothing(prognostic_state(pisces).silicate_climatology)

    run!(simulation)

    picked_up = pisces_simulation(dir)

    picked_up_pisces = picked_up.model.biogeochemistry.underlying_biogeochemistry

    @test all(Array(interior(picked_up_pisces.mixed_layer_depth)) .== 0)

    set!(picked_up; checkpoint = :latest)

    @test picked_up.model.clock.iteration == 2
    @test all(Array(interior(picked_up_pisces.mixed_layer_depth)) .== -25)

    rm(dir; recursive = true)
end

#####
##### Box model
#####

box_model_PAR(t) = 30 * (1 - cos(2π * t / 365days)) + 2

function box_model_simulation(stop_iteration)
    grid = BoxModelGrid()

    clock = Clock(time = zero(grid))

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(FunctionField{Center, Center, Center}(box_model_PAR, grid; clock))

    biogeochemistry = LOBSTER(grid; light_attenuation)

    model = BoxModel(; grid, biogeochemistry, clock)

    set!(model, NO₃ = 10, NH₄ = 0.1, P = 0.1, Z = 0.01)

    simulation = Simulation(model; Δt = 20minutes, stop_iteration, verbose = false)

    simulation.callbacks[:output] = Callback(SpeedyOutput(tempname() * ".jld2"), IterationInterval(1))

    return simulation
end

if architecture == CPU()
    @testset "Checkpointing `BoxModel`" begin
        @info "Testing `BoxModel` checkpointing"

        direct, picked_up = direct_and_picked_up_simulations(box_model_simulation, 5)

        test_matching_models(direct.model, picked_up.model)
    end
end

#####
##### Biogeochemistry without state of its own
#####

@testset "Checkpointing biogeochemistry without state" begin
    grid = RectilinearGrid(architecture; size = (1, 1, 4), extent = (1, 1, 40))

    @test all(isnothing, prognostic_state(LOBSTER(grid)))
    @test all(isnothing, prognostic_state(Biogeochemistry(TracerOnlyBiogeochemistry())))
end

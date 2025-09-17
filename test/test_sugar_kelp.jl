include("dependencies_for_runtests.jl")

using Oceananigans.Architectures: on_architecture

sum_tracer_nitrogen(tracers) = sum(on_architecture(CPU(), interior(tracers.NO₃))) + 
                               sum(on_architecture(CPU(), interior(tracers.NH₄))) +
                               sum(on_architecture(CPU(), interior(tracers.P))) +
                               sum(on_architecture(CPU(), interior(tracers.Z))) +
                               sum(on_architecture(CPU(), interior(tracers.sPON))) +
                               sum(on_architecture(CPU(), interior(tracers.bPON))) +
                               sum(on_architecture(CPU(), interior(tracers.DON)))

sum_tracer_carbon(tracers, redfield, organic_carbon_calcate_ratio) = 
    sum(on_architecture(CPU(), interior(tracers.sPOC))) +
    sum(on_architecture(CPU(), interior(tracers.bPOC))) +
    sum(on_architecture(CPU(), interior(tracers.DOC))) +
    sum(on_architecture(CPU(), interior(tracers.DIC))) + 
    sum(on_architecture(CPU(), interior(tracers.P)) .* (1 + organic_carbon_calcate_ratio) .+ on_architecture(CPU(), interior(tracers.Z))) .* redfield

@testset "SLatissima particle setup and conservations" begin
    grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 1))

    particles = SugarKelpParticles(2; grid, advection = nothing)

    @test particles isa BiogeochemicalParticles
    @test particles.biogeochemistry isa SugarKelp
    @test length(particles) == 2

    biogeochemistry = LOBSTER(grid; 
                              particles, 
                              carbonates = true, 
                              variable_redfield = true, 
                              oxygen = true, 
                              sinking_speeds = NamedTuple())

    model = NonhydrostaticModel(; grid, biogeochemistry, advection = nothing, tracers = (:T, :S))

    set!(model, NO₃ = 10.0, NH₄ = 1.0, DIC = 2000, Alk = 2000, T = 10, S = 35)

    particles.x .= 0.5
    particles.y .= 0.5
    particles.z .= -0.5

    particles.fields.A .= 2
    particles.fields.N .= 1
    particles.fields.C .= 1

    A = on_architecture(CPU(), particles.fields.A)
    N = on_architecture(CPU(), particles.fields.N)
    C = on_architecture(CPU(), particles.fields.C)

    Nₛ = particles.biogeochemistry.structural_nitrogen
    Cₛ = particles.biogeochemistry.structural_carbon
    kₐ = particles.biogeochemistry.structural_dry_weight_per_area

    initial_tracer_N = sum_tracer_nitrogen(model.tracers)
    initial_kelp_N = sum(A .* kₐ .* (N .+ Nₛ)) / (14 * 0.001)

    initial_tracer_C = sum_tracer_carbon(model.tracers, model.biogeochemistry.underlying_biogeochemistry.phytoplankton_redfield, model.biogeochemistry.underlying_biogeochemistry.organic_carbon_calcate_ratio)
    initial_kelp_C = sum(A .* kₐ .* (C .+ Cₛ)) / (12 * 0.001)

    model.clock.time = 60days # get to a high growth phase

    for _ in 1:10
        time_step!(model, 1)
    end

    # not sure we need to repeat this
    A = on_architecture(CPU(), particles.fields.A)
    N = on_architecture(CPU(), particles.fields.N)
    C = on_architecture(CPU(), particles.fields.C)

    final_tracer_N = sum_tracer_nitrogen(model.tracers)
    final_kelp_N = sum(A .* kₐ .* (N .+ Nₛ)) / (14 * 0.001)

    final_tracer_C = sum_tracer_carbon(model.tracers, model.biogeochemistry.underlying_biogeochemistry.phytoplankton_redfield, model.biogeochemistry.underlying_biogeochemistry.organic_carbon_calcate_ratio)
    final_kelp_C = sum(A .* kₐ .* (C .+ Cₛ)) / (12 * 0.001)

    # kelp is being integrated
    @test initial_kelp_N != final_kelp_N
    @test initial_kelp_C != final_kelp_C

    # conservaitons
    # (GPU eps is much larger (~10⁻⁷) than on CPU), is this true??? And is this conservation good enough???
    # this is definitly too high since I didn't catch a sign error of an uptake was the wrong way around
    rtol = ifelse(isa(architecture, CPU), max(√eps(initial_tracer_N + initial_kelp_N), √eps(final_tracer_N + final_kelp_N)), 2e-7)
    @test isapprox(initial_tracer_N + initial_kelp_N, final_tracer_N + final_kelp_N; rtol) 

    rtol = ifelse(isa(architecture, CPU), max(√eps(initial_tracer_C + initial_kelp_C), √eps(final_tracer_C + final_kelp_C)), 7e-7)
    @test isapprox(initial_tracer_C + initial_kelp_C, final_tracer_C + final_kelp_C; rtol)
end
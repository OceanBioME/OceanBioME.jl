using Test, OceanBioME, Oceananigans
using OceanBioME.SLatissimaModel: SLatissima
using Oceananigans.Units
using Oceananigans.Fields: TracerFields
using Oceananigans.Architectures: on_architecture

function intercept_tracer_tendencies!(model, intercepted_tendencies)
    for (name, field) in enumerate(intercepted_tendencies)
        field .= Array(interior(model.timestepper.Gⁿ[name + 3]))
    end
end

sum_tracer_nitrogen(tracers) = sum(Array(interior(tracers.NO₃))) + 
                               sum(Array(interior(tracers.NH₄))) +
                               sum(Array(interior(tracers.P))) +
                               sum(Array(interior(tracers.Z))) +
                               sum(Array(interior(tracers.sPON))) +
                               sum(Array(interior(tracers.bPON))) +
                               sum(Array(interior(tracers.DON)))

sum_tracer_carbon(tracers, redfield, organic_carbon_calcate_ratio) = 
    sum(Array(interior(tracers.sPOC))) +
    sum(Array(interior(tracers.bPOC))) +
    sum(Array(interior(tracers.DOC))) +
    sum(Array(interior(tracers.DIC))) + 
    sum(Array(interior(tracers.P)) .* (1 + organic_carbon_calcate_ratio) .+ Array(interior(tracers.Z))) .* redfield

@testset "SLatissima particle setup and conservations" begin
    grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 1))

    # Initial properties

    particles = SLatissima(; architecture,
                             x = on_architecture(architecture, ones(Float64, 2)),
                             A = on_architecture(architecture, ones(Float64, 2) .* 5),
                             N = on_architecture(architecture, ones(Float64, 2)),
                             C = on_architecture(architecture, ones(Float64, 2)),
                             latitude = 1.0)

    @test length(particles) == 2

    # nitrogen and carbon conservation

    model = NonhydrostaticModel(; grid,
                                  biogeochemistry = LOBSTER(; grid, carbonates = true, 
                                                                    variable_redfield = true, 
                                                                    sinking_speeds = NamedTuple(), 
                                                                    particles),
                                  advection = nothing,
                                  tracers = (:T, :S))

    set!(model, NO₃ = 10.0, NH₄ = 1.0, DIC = 2000, Alk = 2000, T = 10, S = 35)

    initial_tracer_N = sum_tracer_nitrogen(model.tracers)
    initial_kelp_N = sum(Array(particles.A) .* particles.structural_dry_weight_per_area .* (Array(particles.N) .+ particles.structural_nitrogen)) ./ (14 * 0.001)

    initial_tracer_C = sum_tracer_carbon(model.tracers, model.biogeochemistry.underlying_biogeochemistry.phytoplankton_redfield, model.biogeochemistry.underlying_biogeochemistry.organic_carbon_calcate_ratio)
    initial_kelp_C = sum(Array(particles.A) .* particles.structural_dry_weight_per_area .* (Array(particles.C) .+ particles.structural_carbon)) ./ (12 * 0.001)

    model.clock.time = 60days # get to a high growth phase

    for _ in 1:10
        time_step!(model, 1.0)
    end

    final_tracer_N = sum_tracer_nitrogen(model.tracers)
    final_kelp_N = sum(Array(particles.A) .* particles.structural_dry_weight_per_area .* (Array(particles.N) .+ particles.structural_nitrogen)) ./ (14 * 0.001)

    final_tracer_C = sum_tracer_carbon(model.tracers, model.biogeochemistry.underlying_biogeochemistry.phytoplankton_redfield, model.biogeochemistry.underlying_biogeochemistry.organic_carbon_calcate_ratio)
    final_kelp_C = sum(Array(particles.A) .* particles.structural_dry_weight_per_area .* (Array(particles.C) .+ particles.structural_carbon)) ./ (12 * 0.001)

    # kelp is being integrated
    @test initial_kelp_N != final_kelp_N
    @test initial_kelp_C != final_kelp_C

    # conservaitons
    # (GPU eps is much larger (~10⁻⁷) than on CPU)
    rtol = ifelse(isa(architecture, CPU), max(√eps(initial_tracer_N + initial_kelp_N), √eps(final_tracer_N + final_kelp_N)), 2e-7)
    @test isapprox(initial_tracer_N + initial_kelp_N, final_tracer_N + final_kelp_N; rtol) 

    rtol = ifelse(isa(architecture, CPU), max(√eps(initial_tracer_C + initial_kelp_C), √eps(final_tracer_C + final_kelp_C)), 7e-7)
    @test isapprox(initial_tracer_C + initial_kelp_C, final_tracer_C + final_kelp_C; rtol)

    simulation = Simulation(model, Δt = 1.0, stop_iteration = 1)

    # slow but easiest to have this just done on CPU
    intercepted_tendencies = Tuple(Array(interior(field)) for field in values(TracerFields(keys(model.tracers), grid)))

    simulation.callbacks[:intercept_tendencies] = Callback(intercept_tracer_tendencies!; callsite = TendencyCallsite(), parameters = intercepted_tendencies)

    run!(simulation)

    # the model is changing the tracer tendencies - not sure this test actually works as it didn't fail when it should have
    @test any([any(intercepted_tendencies[idx] .!= Array(interior(model.timestepper.Gⁿ[tracer]))) for (idx, tracer) in enumerate((:NO₃, :NH₄, :DIC, :DOC, :bPON, :bPOC))])
end

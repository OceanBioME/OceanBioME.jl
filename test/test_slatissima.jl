using Test, OceanBioME, Oceananigans
using OceanBioME.SLatissimaModel: SLatissima
using Oceananigans.Units
using Oceananigans.Fields: TracerFields

function intercept_tendencies!(model, intercepted_tendencies)
    for tracer in keys(model.tracers)
        intercepted_tendencies[tracer] = Array(interior(model.timestepper.Gⁿ[tracer]))
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

    particles = SLatissima(; x = ones(Float64, 2),
                             A = ones(Float64, 2) .* 5,
                             N = ones(Float64, 2),
                             C = ones(Float64, 2),
                             latitude = 1.0,
                             pescribed_temperature = (args...) -> 10.0,
                             pescribed_salinity = (args...) -> 35.0)

    @test length(particles) == 2

    # nitrogen and carbon conservation

    model = NonhydrostaticModel(; grid,
                                  biogeochemistry = LOBSTER(; grid, carbonates = true, 
                                                                    variable_redfield = true, 
                                                                    sinking_speeds = NamedTuple(), 
                                                                    particles),
                                  advection = nothing)

    set!(model, NO₃ = 10.0, NH₄ = 1.0, DIC = 2000, Alk = 2000)

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
    @test initial_tracer_N + initial_kelp_N ≈ final_tracer_N + final_kelp_N
    @test initial_tracer_C + initial_kelp_C ≈ final_tracer_C + final_kelp_C

    simulation = Simulation(model, Δt = 1.0, stop_iteration = 1)

    # slow but easiest to have this just done on CPU
    intercepted_tendencies = Tuple(Array(interior(field)) for field in values(TracerFields(keys(model.tracers), grid)))

    simulation.callbacks[:intercept_tendencies] = Callback(intercept_tendencies!; callsite = TendencyCallsite(), parameters = intercepted_tendencies)

    run!(simulation)

    # the model is changing the tracer tendencies - not sure this test actually works as it didn't fail when it should have
    @test any([any(intercepted_tendencies[idx] .!= Array(interior(model.timestepper.Gⁿ[tracer]))) for (idx, tracer) in enumerate((:NO₃, :NH₄, :DIC, :DOC, :bPON, :bPOC))])
end

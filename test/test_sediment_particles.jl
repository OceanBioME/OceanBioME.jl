using OceanBioME, Oceananigans, Test, JLD2

using OceanBioME.Sediments: SimpleMultiG
using Oceananigans.Units

using OceanBioME.Sediments: sediment_tracers, sediment_fields
using Oceananigans: Field
using Oceananigans.Fields: TracerFields

function intercept_tendencies!(model, intercepted_tendencies)
    for tracer in keys(model.tracers)
        copyto!(intercepted_tendencies[tracer], model.timestepper.Gⁿ[tracer])
    end
end

function test_flat_sediment(architecture; timestepper = :QuasiAdamsBashforth2)
    grid = RectilinearGrid(architecture; size=(3, 3, 3), extent=(10, 10, 200))

    sediment_model = SimpleMultiG(; grid)

    particles = SLatissima(; x = ones(Float64, 2), 
                             A = ones(Float64, 2) .* 5, 
                             N = ones(Float64, 2), 
                             C = ones(Float64, 2), 
                             latitude = 1.0,
                             pescribed_temperature = (args...) -> 10.0,
                             pescribed_salinity = (args...) -> 35.0)

    biogeochemistry = LOBSTER(; grid, 
                                carbonates = true, oxygen = true, variable_redfield = true, 
                                open_bottom = true, 
                                surface_phytosynthetically_active_radiation = (x, y, t) -> 80,
                                sediment_model,
                                particles)

    model = NonhydrostaticModel(;grid, biogeochemistry, 
                                 closure = nothing,
                                 timestepper)

    set!(model.biogeochemistry.sediment_model.fields.N_fast, 0.0230)
    set!(model.biogeochemistry.sediment_model.fields.N_slow, 0.0807)

    set!(model.biogeochemistry.sediment_model.fields.C_fast, 0.5893)
    set!(model.biogeochemistry.sediment_model.fields.C_slow, 0.1677)

    set!(model, P = 0.4686, Z = 0.5363, 
                NO₃ = 2.3103, NH₄ = 0.0010, 
                DIC = 2106.9, Alk = 2408.9, 
                O₂ = 258.92, 
                DOC = 5.3390, DON = 0.8115,
                sPON = 0.2299, sPOC = 1.5080,
                bPON = 0.0103, bPOC = 0.0781)

    initial_kelp_N = sum(particles.A .* particles.structural_dry_weight_per_area .* (particles.N .+ particles.structural_nitrogen)) ./ (14 * 0.001)

    simulation = Simulation(model, Δt = 1, stop_iteration = 1)

    intercepted_tendencies = TracerFields(keys(model.tracers), grid)

    simulation.callbacks[:intercept_tendencies] = Callback(intercept_tendencies!; callsite = TendencyCallsite(), parameters = intercepted_tendencies)

    run!(simulation)

    final_kelp_N = sum(particles.A .* particles.structural_dry_weight_per_area .* (particles.N .+ particles.structural_nitrogen)) ./ (14 * 0.001)

    # kelp is getting integrated
    @test initial_kelp_N != final_kelp_N

    # the model is changing the tracer tendencies
    @test any([any(intercepted_tendencies[tracer] .!= model.timestepper.Gⁿ[tracer]) for tracer in keys(model.tracers)])

    # the sediment tendencies are being updated
    @test all([any(tend .!= 0.0) for tend in model.biogeochemistry.sediment_model.tendencies.Gⁿ])
    @test all([any(tend .!= 0.0) for tend in model.biogeochemistry.sediment_model.tendencies.G⁻])

    # the sediment values are being integrated
    initial_values = (N_fast = 0.0230, N_slow = 0.0807, C_fast = 0.5893, C_slow = 0.1677, N_ref = 0.0, C_ref = 0.0)
    @test all([any(field .!= initial_values[name]) for (name, field) in pairs(model.biogeochemistry.sediment_model.fields)])

    # kelp is being integrated
    @test initial_kelp_N != final_kelp_N
    return nothing
end

@testset "Sediment" begin
    for arch in (CPU(), )
        for timestepper in (:QuasiAdamsBashforth2, :RungeKutta3)
            @info "Testing sediment on $arch with $timestepper"
            @testset "$arch, $timestepper" test_flat_sediment(arch; timestepper)
        end
    end
end

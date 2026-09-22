include("dependencies_for_runtests.jl")

# using OceanBioME: setup_velocity_fields, valid_sinking_velocity_locations

# using Oceananigans.Fields: AbstractField, CenterField, ConstantField, FunctionField, ZFaceField, location

function with_negative_tracers(bgc, negative_tracers)
    return Biogeochemistry(bgc.underlying_biogeochemistry;
                           light_attenuation = bgc.light_attenuation,
                           sediment = bgc.sediment,
                           particles = bgc.particles,
                           modifiers = bgc.modifiers,
                           negative_tracers)
end


struct ContinuousNegativeTracerTestBGC <: Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry end

Oceananigans.Biogeochemistry.required_biogeochemical_tracers(::ContinuousNegativeTracerTestBGC) = (:A, :T)

@inline (::ContinuousNegativeTracerTestBGC)(::Val{:A}, x, y, z, t, A, T, auxiliary) = A + T + auxiliary

struct DiscreteNegativeTracerTestBGC <: Oceananigans.Biogeochemistry.AbstractBiogeochemistry end

Oceananigans.Biogeochemistry.required_biogeochemical_tracers(::DiscreteNegativeTracerTestBGC) = (:A, :T)
Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(::DiscreteNegativeTracerTestBGC) = NamedTuple()

@inline function (::DiscreteNegativeTracerTestBGC)(i, j, k, grid, ::Val{:A}, clock, fields, auxiliary_fields)
    return fields[:A][i, j, k] + fields.T[i, j, k] + fields.extra[i, j, k]
end

function test_negative_scaling(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    bgc = NPZD(grid)
    scaler = ScaleNegativeTracers(bgc.underlying_biogeochemistry)
    model = NonhydrostaticModel(grid; biogeochemistry = with_negative_tracers(bgc, scaler))

    set!(model, N = 2, P = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)
    run!(simulation)

    N = Array(interior(model.tracers.N))[1, 1, 1]
    P = Array(interior(model.tracers.P))[1, 1, 1]

    return (N ≈ 1) && (P ≈ 0.0)
end

function test_negative_clipping(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    bgc = NPZD(grid)
    clip = ClipNegativeTracers(; exclude = (:Z, ))
    model = NonhydrostaticModel(grid; biogeochemistry = with_negative_tracers(bgc, clip))

    set!(model, N = 2, P = -1, Z = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)
    run!(simulation)

    N = Array(interior(model.tracers.N))[1, 1, 1]
    P = Array(interior(model.tracers.P))[1, 1, 1]
    Z = Array(interior(model.tracers.Z))[1, 1, 1]

    return (N ≈ 2) && (P ≈ 0.0) && (Z ≈ -1)
end

@testset "Test negative tracer handling" begin
    @test test_negative_scaling(architecture)
    @test test_negative_clipping(architecture)

    grid = RectilinearGrid(architecture, size = (1, 1, 1), extent = (1, 1, 1))
    bgc = NPZD(grid)
    @test isnothing(bgc.negative_tracers)

    raw = Biogeochemistry(ContinuousNegativeTracerTestBGC())
    safe = Biogeochemistry(ContinuousNegativeTracerTestBGC(); negative_tracers = IgnoreNegativeTracerValues())

    @test raw(Val(:A), 0, 0, 0, 0, -1.0, -2.0, -3.0) == -6.0
    @test safe(Val(:A), 0, 0, 0, 0, -1.0, -2.0, -3.0) == -5.0
    @test safe(Val(:A), 0, 0, 0, 0, 1.0, -2.0, -3.0) == raw(Val(:A), 0, 0, 0, 0, 1.0, -2.0, -3.0)

    raw_discrete = Biogeochemistry(DiscreteNegativeTracerTestBGC())
    safe_discrete = Biogeochemistry(DiscreteNegativeTracerTestBGC(); negative_tracers = IgnoreNegativeTracerValues())
    fields = (A = fill(-1.0, 1, 1, 1), T = fill(-2.0, 1, 1, 1), extra = fill(-3.0, 1, 1, 1))

    @test raw_discrete(1, 1, 1, grid, Val(:A), nothing, fields) == -6.0
    @test safe_discrete(1, 1, 1, grid, Val(:A), nothing, fields) == -5.0
    @test fields.A[1, 1, 1] == -1.0

    light = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(10.0))
    reference = NPZD(grid; light_attenuation = light)
    safe_npzd = with_negative_tracers(NPZD(grid; light_attenuation = light), IgnoreNegativeTracerValues())
    reference_model = NonhydrostaticModel(grid; biogeochemistry = reference)
    safe_model = NonhydrostaticModel(grid; biogeochemistry = safe_npzd)

    common = (; N = 2.0, Z = 0.1, D = 0.1, T = -1.0)
    set!(reference_model; common..., P = 0.0)
    set!(safe_model; common..., P = -1e-3)

    tendency_reference = reference(1, 1, 1, grid, Val(:P), reference_model.clock, reference_model.tracers)
    tendency_safe = safe_npzd(1, 1, 1, grid, Val(:P), safe_model.clock, safe_model.tracers)

    @test tendency_safe ≈ tendency_reference
    @test safe_model.tracers.P[1, 1, 1] == -1e-3

    for light_model in (TwoBandPhotosyntheticallyActiveRadiation,
                        MultiBandPhotosyntheticallyActiveRadiation)
        reference = NPZD(grid; light_attenuation = light_model(grid, 100.0))
        safe_npzd = with_negative_tracers(NPZD(grid; light_attenuation = light_model(grid, 100.0)),
                                          IgnoreNegativeTracerValues())
        reference_model = NonhydrostaticModel(grid; biogeochemistry = reference)
        safe_model = NonhydrostaticModel(grid; biogeochemistry = safe_npzd)

        set!(reference_model; common..., P = 0.0)
        set!(safe_model; common..., P = -1e-3)

        @test OceanBioME.chlorophyll(safe_npzd, safe_model)[1, 1, 1] == 0

        Oceananigans.Biogeochemistry.update_biogeochemical_state!(reference_model, reference.light_attenuation)
        Oceananigans.Biogeochemistry.update_biogeochemical_state!(safe_model, safe_npzd.light_attenuation)

        reference_light = Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(reference.light_attenuation)
        safe_light = Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(safe_npzd.light_attenuation)

        @test all(name -> safe_light[name][1, 1, 1] ≈ reference_light[name][1, 1, 1], keys(reference_light))
        @test all(name -> isfinite(safe_light[name][1, 1, 1]), keys(safe_light))
        @test safe_model.tracers.P[1, 1, 1] == -1e-3
    end
end

# scalar_sinking_speeds = (A = 1, B = 1.0)
#
# grid = RectilinearGrid(architecture, size = (1, 1, 10), extent = (1, 1, 10))
#
# field_sinking_speeds = (C = ConstantField(1), D = FunctionField{Center, Center, Face}((x, y, z) -> z, grid), E = ZFaceField(grid))
#
# sinking_speeds = merge(scalar_sinking_speeds, field_sinking_speeds)
#
# @testset "Test sinking velocity setup" begin
#     sinking_velocities = @test_nowarn setup_velocity_fields(sinking_speeds, grid, true)
#
#     @test all(map(w -> isa(w, AbstractField) & (location(w) in valid_sinking_velocity_locations), values(sinking_velocities)))
#
#     sinking_velocities =
#         @test_warn ("The sinking velocity provided for C is a field and therefore `open_bottom=false` can't be enforced automatically",
#                     "The sinking velocity provided for D is a field and therefore `open_bottom=false` can't be enforced automatically",
#                     "The sinking velocity provided for E is a field and therefore `open_bottom=false` can't be enforced automatically") setup_velocity_fields(sinking_speeds, grid, false)
#
#     @test all([isa(w, AbstractField) & (location(w) in valid_sinking_velocity_locations) for w in values(sinking_velocities)])
#
#     @test all(map(w -> Array(interior(w, 1, 1, 1)) .== 0, sinking_velocities[(:A, :B)]))
#
#     @test_warn "The location of the sinking velocity field provided for X is incorrect, it should be (Center, Center, Face)" setup_velocity_fields((X = CenterField(grid), ), grid, true)
# end

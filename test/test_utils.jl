include("dependencies_for_runtests.jl")

using OceanBioME: setup_velocity_fields, valid_sinking_velocity_locations

using Oceananigans.Architectures: on_architecture
using Oceananigans.AbstractOperations: AbstractOperation
using Oceananigans.Fields: AbstractField, CenterField, ConstantField, FunctionField, ZFaceField, location

scalar_value(field::AbstractField) = on_architecture(CPU(), interior(field, 1, 1, 1))[1]

function scalar_value(operation::AbstractOperation)
    field = Field(operation)
    compute!(field)
    return scalar_value(field)
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

    bgc = NPZD(grid; negative_tracers = ScaleNegativeTracers())
    model = NonhydrostaticModel(grid; biogeochemistry = bgc)

    set!(model, N = 2, P = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)
    run!(simulation)

    N = scalar_value(model.tracers.N)
    P = scalar_value(model.tracers.P)

    return (N ≈ 1) && (P ≈ 0.0)
end

function test_negative_clipping(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    bgc = NPZD(grid; negative_tracers = ClipNegativeTracers(; exclude = (:Z, )))
    model = NonhydrostaticModel(grid; biogeochemistry = bgc)

    set!(model, N = 2, P = -1, Z = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)
    run!(simulation)

    N = scalar_value(model.tracers.N)
    P = scalar_value(model.tracers.P)
    Z = scalar_value(model.tracers.Z)

    return (N ≈ 2) && (P ≈ 0.0) && (Z ≈ -1)
end

@testset "Test negative tracer handling" begin
    @test test_negative_scaling(architecture)
    @test test_negative_clipping(architecture)

    grid = RectilinearGrid(architecture, size = (1, 1, 1), extent = (1, 1, 1))
    bgc = NPZD(grid)
    ignore = NPZD(grid; negative_tracers = IgnoreNegativeTracerValues()).negative_tracers
    scale = NPZD(grid; negative_tracers = ScaleNegativeTracers()).negative_tracers
    scale_compat = NPZD(grid; scale_negatives = true).negative_tracers
    composed_npzd = NPZD(grid; negative_tracers = (IgnoreNegativeTracerValues(), ScaleNegativeTracers()))

    @test isnothing(bgc.negative_tracers)
    @test ignore isa IgnoreNegativeTracerValues
    @test scale isa ScaleNegativeTracers
    @test scale_compat isa ScaleNegativeTracers
    @test composed_npzd.negative_tracers[1] isa IgnoreNegativeTracerValues
    @test composed_npzd.negative_tracers[2] isa ScaleNegativeTracers
    @test_throws ArgumentError NPZD(grid; negative_tracers = IgnoreNegativeTracerValues(), scale_negatives = true)

    raw = Biogeochemistry(ContinuousNegativeTracerTestBGC())
    safe = Biogeochemistry(ContinuousNegativeTracerTestBGC(); negative_tracers = IgnoreNegativeTracerValues())

    @test raw(Val(:A), 0, 0, 0, 0, -1.0, -2.0, -3.0) == -6.0
    @test safe(Val(:A), 0, 0, 0, 0, -1.0, -2.0, -3.0) == -5.0
    @test safe(Val(:A), 0, 0, 0, 0, 1.0, -2.0, -3.0) == raw(Val(:A), 0, 0, 0, 0, 1.0, -2.0, -3.0)

    composed = Biogeochemistry(ContinuousNegativeTracerTestBGC();
                                negative_tracers = (ClipNegativeTracers(), (IgnoreNegativeTracerValues(),)))
    @test composed(Val(:A), 0, 0, 0, 0, -1.0, -2.0, -3.0) == -5.0

    raw_discrete = Biogeochemistry(DiscreteNegativeTracerTestBGC())
    safe_discrete = Biogeochemistry(DiscreteNegativeTracerTestBGC(); negative_tracers = IgnoreNegativeTracerValues())
    fields = (A = fill(-1.0, 1, 1, 1), T = fill(-2.0, 1, 1, 1), extra = fill(-3.0, 1, 1, 1))

    @test raw_discrete(1, 1, 1, grid, Val(:A), nothing, fields) == -6.0
    @test safe_discrete(1, 1, 1, grid, Val(:A), nothing, fields) == -5.0

    composed_discrete = Biogeochemistry(DiscreteNegativeTracerTestBGC();
                                         negative_tracers = (ClipNegativeTracers(), IgnoreNegativeTracerValues()))
    @test composed_discrete(1, 1, 1, grid, Val(:A), nothing, fields) == -5.0
    @test fields.A[1, 1, 1] == -1.0

    light = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(10.0))
    reference = NPZD(grid; light_attenuation = light)
    safe_npzd = NPZD(grid; light_attenuation = light, negative_tracers = IgnoreNegativeTracerValues())
    reference_model = NonhydrostaticModel(grid; biogeochemistry = reference)
    safe_model = NonhydrostaticModel(grid; biogeochemistry = safe_npzd)

    common = (; N = 2.0, Z = 0.1, D = 0.1, T = -1.0)
    set!(reference_model; common..., P = 0.0)
    set!(safe_model; common..., P = -1e-3)

    tendency_reference = CUDA.@allowscalar reference(1, 1, 1, grid, Val(:P), reference_model.clock, reference_model.tracers)
    tendency_safe = CUDA.@allowscalar safe_npzd(1, 1, 1, grid, Val(:P), safe_model.clock, safe_model.tracers)

    @test tendency_safe ≈ tendency_reference

    simulation = Simulation(safe_model, Δt = 1e-10, stop_iteration = 1)
    run!(simulation)
    @test scalar_value(safe_model.tracers.P) ≈ -1e-3

    for light_model in (TwoBandPhotosyntheticallyActiveRadiation,
                        MultiBandPhotosyntheticallyActiveRadiation)
        reference = NPZD(grid; light_attenuation = light_model(grid, 100.0))
        safe_npzd = NPZD(grid;
                         light_attenuation = light_model(grid, 100.0),
                         negative_tracers = IgnoreNegativeTracerValues())
        reference_model = NonhydrostaticModel(grid; biogeochemistry = reference)
        safe_model = NonhydrostaticModel(grid; biogeochemistry = safe_npzd)

        set!(reference_model; common..., P = 0.0)
        set!(safe_model; common..., P = -1e-3)


        Oceananigans.Biogeochemistry.update_biogeochemical_state!(reference_model, reference.light_attenuation)
        Oceananigans.Biogeochemistry.update_biogeochemical_state!(safe_model, safe_npzd.light_attenuation)

        reference_light = Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(reference.light_attenuation)
        safe_light = Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(safe_npzd.light_attenuation)

        @test all(name -> scalar_value(safe_light[name]) ≈ scalar_value(reference_light[name]), keys(reference_light))
        @test all(name -> isfinite(scalar_value(safe_light[name])), keys(safe_light))
        @test scalar_value(safe_model.tracers.P) == -1e-3
    end
end

scalar_sinking_speeds = (A = 1, B = 1.0)

grid = RectilinearGrid(architecture, size = (1, 1, 10), extent = (1, 1, 10))

field_sinking_speeds = (C = ConstantField(1), D = FunctionField{Center, Center, Face}((x, y, z) -> z, grid), E = ZFaceField(grid))

sinking_speeds = merge(scalar_sinking_speeds, field_sinking_speeds)

@testset "Test sinking velocity setup" begin
    sinking_velocities = @test_nowarn setup_velocity_fields(sinking_speeds, grid, true)

    @test all(map(w -> isa(w, AbstractField) & (location(w) in valid_sinking_velocity_locations), values(sinking_velocities)))

    sinking_velocities =
        @test_warn ("The sinking velocity provided for C is a field and therefore `open_bottom=false` can't be enforced automatically",
                    "The sinking velocity provided for D is a field and therefore `open_bottom=false` can't be enforced automatically",
                    "The sinking velocity provided for E is a field and therefore `open_bottom=false` can't be enforced automatically") setup_velocity_fields(sinking_speeds, grid, false)

    @test all([isa(w, AbstractField) & (location(w) in valid_sinking_velocity_locations) for w in values(sinking_velocities)])

    @test all(map(w -> Array(interior(w, 1, 1, 1)) .== 0, sinking_velocities[(:A, :B)]))

    @test_warn "The location of the sinking velocity field provided for X is incorrect, it should be (Center, Center, Face)" setup_velocity_fields((X = CenterField(grid), ), grid, true)
end

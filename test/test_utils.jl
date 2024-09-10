include("dependencies_for_runtests.jl")

using OceanBioME: setup_velocity_fields, valid_sinking_velocity_locations

using Oceananigans.Fields: AbstractField, CenterField, ConstantField, FunctionField, ZFaceField

function test_negative_scaling(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(; grid, biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid, scale_negatives = true))

    set!(model, N = 2, P = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)

    run!(simulation)

    N = Array(interior(model.tracers.N))[1, 1, 1]
    P = Array(interior(model.tracers.P))[1, 1, 1]

    return (N ≈ 1) && (P ≈ 0.0)
end

function test_negative_zeroing(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(; grid, biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid, modifiers = ZeroNegativeTracers(; exclude = (:Z, ))))

    set!(model, N = 2, P = -1, Z = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)

    run!(simulation)

    N = Array(interior(model.tracers.N))[1, 1, 1]
    P = Array(interior(model.tracers.P))[1, 1, 1]
    Z = Array(interior(model.tracers.Z))[1, 1, 1]

    return (N ≈ 2) && (P ≈ 0.0) && (Z ≈ -1)
end

@testset "Test negative tracer handeling" begin
    @test test_negative_scaling(architecture)
    @test test_negative_zeroing(architecture)
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

    @test all(map(w -> Array(interior(w, 1, 1, 1)) == 0, sinking_velocities[(:A, :B)]))

    @test_warn "The location of the sinking velocity field provided for X is incorrect, it should be (Center, Center, Face)" setup_velocity_fields((X = CenterField(grid), ), grid, true)
end

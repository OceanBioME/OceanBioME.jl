include("dependencies_for_runtests.jl")

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

@testset "Test Utils" begin
    @test test_negative_scaling(architecture)
    @test test_negative_zeroing(architecture)
end

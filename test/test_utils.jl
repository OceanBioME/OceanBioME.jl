using Test, OceanBioME, Oceananigans

function test_column_diffusion_timescale(arch)
    κ = 1e-3
    @inline κₜ(x, y, z, t) = κ

    grid = RectilinearGrid(arch, size=(1, 1, 5), x=(0, 10), y=(0, 3), z=[-1, -0.6, -0.5, -0.2, -0.19, 0])
    min_Δz = minimum_zspacing(grid)

    model = NonhydrostaticModel(; grid,
                                closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
                                biogeochemistry = LOBSTER(; grid,
                                                            surface_phytosynthetically_active_radiation = (x, y, t) -> 100,
                                                            carbonates = true),
                                advection = nothing)

    return column_diffusion_timescale(model) ≈ min_Δz ^ 2 / κ
end

function test_negative_scaling(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(; grid, biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid, scale_negatives = true))

    set!(model, N = 2, P = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)

    run!(simulation)

    return (model.tracers.N[1, 1, 1] ≈ 1) && (model.tracers.P[1, 1, 1] ≈ 0.0)
end

function test_negative_zeroing(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(; grid, biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid, modifiers = ZeroNegativeTracers(; exclude = (:Z, ))))

    set!(model, N = 2, P = -1, Z = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)

    run!(simulation)

    return (model.tracers.N[1, 1, 1] ≈ 2) && (model.tracers.P[1, 1, 1] ≈ 0.0) && (model.tracers.Z[1, 1, 1] ≈ -1)
end

@testset "Test Utils" begin
    @test test_column_diffusion_timescale(architecture)
    @test test_negative_scaling(architecture)
    @test test_negative_zeroing(architecture)
end

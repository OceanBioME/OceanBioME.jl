using Test, OceanBioME, Oceananigans

function test_column_diffusion_timescale(arch)
    κ = 1e-3
    @inline κₜ(x, y, z, t) = κ

    grid = RectilinearGrid(arch, size=(1, 1, 5), x=(0, 10), y=(0, 3), z=[-1, -0.6, -0.5, -0.2, -0.19, 0])
    min_Δz = 1e-2 # for the grid above

    model = NonhydrostaticModel(; grid,
                                closure = ScalarDiffusivity(ν = κₜ, κ = κₜ), 
                                biogeochemistry = LOBSTER(; grid,
                                                            surface_phytosynthetically_active_radiation = (x, y, t) -> 100,
                                                            carbonates = true),
                                advection = nothing)

    return column_diffusion_timescale(model) ≈ min_Δz^2 / κ
end

function test_negative_scaling(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(; grid, tracers = (:A, :B))

    set!(model, A = 1, B = -0.5)

    simulation = Simulation(model, Δt = 1)

    negativity_protection! = ScaleNegativeTracers(; model, tracers = (:A, :B))
    
    negativity_protection!(simulation)

    return (model.tracers.A[1, 1, 1] == 0.5) && (model.tracers.B[1, 1, 1] == 0.0)
end

function test_negative_scaling(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(; grid, tracers = :A)

    model.timestepper.Gⁿ.A .= NaN

    remove_NaN_tendencies!(model)

    return model.timestepper.Gⁿ.A[1, 1, 1] == 0.0
end

function test_negative_zeroing(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(; grid, tracers = (:A, :B))

    set!(model, A = -1, B = -0.5)

    zero_negative_tracers!(model, params = (exclude = (:B, ), ))

    return (model.tracers.A[1, 1, 1] == 0.0) && (model.tracers.B[1, 1, 1] == -0.5)
end

@testset "Test Utils" begin
    for arch in (CPU(), )
        @test test_column_diffusion_timescale(arch)
        @test test_negative_scaling(arch)
        @test test_negative_scaling(arch)
        @test test_negative_zeroing(arch)
    end
end

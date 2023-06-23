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

@testset "Test Utils" begin
    for arch in (CPU(), )
        @test test_column_diffusion_timescale(arch)
    end
end

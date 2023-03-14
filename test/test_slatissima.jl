using Test
using OceanBioME: SLatissima
using CUDA

@testset "Sugar Kelp Particle Setup" begin
    # Initial properties
    kelp_particles = SLatissima.setup(n = 1, x₀ = 1.0, y₀ = 1.0, z₀ = 1.0, A₀ = 1.0, N₀ = 1.0, C₀ = 1.0, latitude = 1.0)

    @test length(kelp_particles) == 1
    @test all(kelp_particles.properties.x .== 1.0)

    kelp_particles = SLatissima.setup(n = 2, x₀ = 1.0, y₀ = 1.0, z₀ = 1.0, A₀ = 1.0, N₀ = 1.0, C₀ = 1.0, latitude = 1.0)

    @test length(kelp_particles) == 2
    @test all(kelp_particles.properties.x .== 1.0)

    kelp_particles = SLatissima.setup(n = 2, x₀ = [1.0, 2.0], y₀ = 1.0, z₀ = 1.0, A₀ = 1.0, N₀ = 1.0, C₀ = 1.0, latitude = 1.0)

    @test length(kelp_particles) == 2
    @test all(kelp_particles.properties.x == [1.0, 2.0])
    @test all(kelp_particles.properties.y .== 1.0)

    # Temperature and salinity function passing

    @test !(:T in keys(kelp_particles.parameters.equation_parameters)) && !(:S in keys(kelp_particles.parameters.equation_parameters))

    kelp_particles = SLatissima.setup(n = 1, x₀ = 1.0, y₀ = 1.0, z₀ = 1.0, A₀ = 1.0, N₀ = 1.0, C₀ = 1.0, latitude = 1.0, T = (x, y, z, t) -> 1.0, S = (x, y, z, t) -> 35.0, urel = 0.2)
    @test (:T in keys(kelp_particles.parameters.equation_parameters)) && (:S in keys(kelp_particles.parameters.equation_parameters))

    # Optional tracers

    @test kelp_particles.parameters.coupled_fields == (NO₃ = :j_NO₃,)
    @test kelp_particles.parameters.tracked_fields == (NO₃ = :NO₃, PAR = :PAR)

    kelp_particles = SLatissima.setup(n = 1, x₀ = 1.0, y₀ = 1.0, z₀ = 1.0, A₀ = 1.0, N₀ = 1.0, C₀ = 1.0, latitude = 1.0; optional_tracers = (:NH₄, :DIC, :bPON, :bPOC, :O₂, :DON, :DOC))

    @test kelp_particles.parameters.tracked_fields == (NO₃ = :NO₃, PAR = :PAR, NH₄ = :NH₄, T = :T, S = :S, u = :u, v = :v, w = :w)
    @test kelp_particles.parameters.coupled_fields == (NO₃ = :j_NO₃, NH₄ = :j_NH₄, DIC = :j_DIC, bPON = :νⁿ, bPOC = :νᶜ, O₂ = :j_OXY, DON = :eⁿ, DOC = :eᶜ)

    return nothing
end

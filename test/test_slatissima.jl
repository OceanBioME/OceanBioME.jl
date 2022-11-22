using Test
using OceanBioME: SLatissima

@testset "Sugar Kelp Particle Setup" begin
    # Initial properties

    kelp_particles = SLatissima.setup(1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    @test length(kelp_particles) == 1
    @test all(kelp_particles.properties.x .== 1.0)

    kelp_particles = SLatissima.setup(2, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    @test length(kelp_particles) == 2
    @test all(kelp_particles.properties.x .== 1.0)

    kelp_particles = SLatissima.setup(2, [1.0, 2.0], 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

    @test length(kelp_particles) == 2
    @test all(kelp_particles.properties.x == [1.0, 2.0])
    @test all(kelp_particles.properties.y .== 1.0)

    # Temperature and salinity function passing

    @test !(:T in keys(kelp_particles.parameters.equation_parameters)) && !(:S in keys(kelp_particles.parameters.equation_parameters))

    kelp_particles = SLatissima.setup(1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0; T = (x, y, z, t) -> 1.0, S = (x, y, z, t) -> 35.0, urel = 0.2)
    @test (:T in keys(kelp_particles.parameters.equation_parameters)) && (:S in keys(kelp_particles.parameters.equation_parameters))

    # Optional tracers

    @test kelp_particles.parameters.coupled_fields == (NO₃ = :j_NO₃,)
    @test kelp_particles.parameters.tracked_fields == (NO₃ = :NO₃, PAR = :PAR)

    kelp_particles = SLatissima.setup(1, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0; optional_tracers = (:NH₄, :DIC, :DD, :DDᶜ, :OXY, :DOM))

    @test kelp_particles.parameters.tracked_fields == (NO₃ = :NO₃, PAR = :PAR, NH₄ = :NH₄)
    @test kelp_particles.parameters.coupled_fields == (NO₃ = :j_NO₃, NH₄ = :j_NH₄, DIC = :j_DIC, DD = :νⁿ, DDᶜ = :νᶜ, OXY = :j_OXY, DOM = :e)
end
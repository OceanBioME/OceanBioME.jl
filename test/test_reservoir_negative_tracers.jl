using OceanBioME: conserved_tracers
using OceanBioME: reservoir_negative_values
using Oceananigans.TimeSteppers: update_state!

@testset "ReservoirNegativeTracers" begin
    # Pure-cell calculation: a negative total survives as a signed reservoir,
    # and the correction is idempotent when there is nothing left to repay it.
    for FT in (Float32, Float64)
        weights = ((one(FT), one(FT)),)
        corrected, reservoirs = reservoir_negative_values((FT(0.25), FT(-1)),
                                                          (zero(FT),), weights, (1,))
        @test corrected == (zero(FT), zero(FT))
        @test reservoirs[1] ≈ FT(-0.75)
        @test sum((FT(0.25), FT(-1))) ≈ sum(corrected) + reservoirs[1]
        @test reservoir_negative_values(corrected, reservoirs, weights, (1,)) == (corrected, reservoirs)
    end

    grid = RectilinearGrid(architecture, size = (1, 1, 1), extent = (1, 1, 1))
    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))
    inorganic_carbon = CarbonateSystem()
    base = NPZD(grid; inorganic_carbon, light_attenuation)
    groups = conserved_tracers(base)
    modifier = ReservoirNegativeTracers(base; repayment_tracers = (nitrogen = :N, carbon = :DIC))
    biogeochemistry = NPZD(grid; inorganic_carbon, light_attenuation, modifiers = modifier)
    model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing)

    @test (:nitrogen_reservoir in keys(model.tracers)) && (:carbon_reservoir in keys(model.tracers))
    @test OceanBioME.is_reservoir_tendency(modifier, 1, 1, 1, grid, Val(:nitrogen_reservoir))
    @test !OceanBioME.is_reservoir_tendency(modifier, 1, 1, 1, grid, Val(:N))

    for values in ((N = 2.0, P = -0.5, Z = 1.0, D = 1.0, DIC = 100.0),
                   (N = 0.1, P = -1.0, Z = 0.0, D = 0.0, DIC = 0.0))
        set_tracers!(model; N = 0, P = 0, Z = 0, D = 0, DIC = 0,
                            nitrogen_reservoir = 0, carbon_reservoir = 0)
        set_tracers!(model; values...)

        before = group_totals(model, groups)
        update_state!(model)
        after = group_totals(model, groups)
        after = (nitrogen = after.nitrogen + tracer_value(model, :nitrogen_reservoir),
                carbon = after.carbon + tracer_value(model, :carbon_reservoir))

        @test all(map(≈, before, after))
        @test all(name -> tracer_value(model, name) >= 0, (:N, :P, :Z, :D, :DIC))
        @test tracer_value(model, :nitrogen_reservoir) <= 0
        @test tracer_value(model, :carbon_reservoir) <= 0
    end

    # The reservoir can repay a later positive inorganic pool without changing
    # the two independently conserved element inventories.
    set_tracers!(model; N = 2, P = 0, Z = 0, D = 0, DIC = 10,
                        nitrogen_reservoir = -0.5, carbon_reservoir = -1)
    update_state!(model)
    @test tracer_value(model, :N) ≈ 1.5
    @test tracer_value(model, :DIC) ≈ 9
    @test tracer_value(model, :nitrogen_reservoir) ≈ 0
    @test tracer_value(model, :carbon_reservoir) ≈ 0

    @test_throws ArgumentError ReservoirNegativeTracers(
        (nitrogen = (N = 1, P = 1), carbon = (DIC = 1, P = 6.56));
        repayment_tracers = (nitrogen = :P, carbon = :DIC))
    @test_throws ArgumentError ReservoirNegativeTracers(
        (oxygen = (O₂ = 1, P = -7),); repayment_tracers = (oxygen = :O₂,))
end

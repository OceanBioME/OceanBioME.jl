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

@testset "Reservoir diffusion and repayment" begin
    # An abiotic, single-nutrient model isolates Oceananigans transport from
    # biological reactions. The reservoir is an ordinary prognostic tracer.
    grid = RectilinearGrid(architecture; size = 16, x = (0, 16),
                           topology = (Bounded, Flat, Flat))
    base = NPZD(grid; plankton = Abiotic(), detritus = InstantRemineralisationDetritus(),
                light_attenuation = nothing)
    modifier = ReservoirNegativeTracers(base.underlying_biogeochemistry;
                                       repayment_tracers = (nitrogen = :N,))
    bgc = Biogeochemistry(base.underlying_biogeochemistry; modifiers = modifier)

    closure = HorizontalScalarDiffusivity(ν = 0, κ = (N = 0, nitrogen_reservoir = 0.5))
    model = NonhydrostaticModel(grid; biogeochemistry = bgc, closure,
                                momentum_advection = nothing,
                                tracer_advection = (N = Centered(), nitrogen_reservoir = nothing))

    # Cell 1 is negative, cell 2 contains no nutrient, and cells 3+ can
    # repay the diffusing debt. Boundaries are no-flux.
    set!(model.tracers.N, (x, y, z) -> x < 1 ? -0.5 : x < 2 ? 0.0 : 1.0)
    update_state!(model)

    concentrations(field) = vec(Array(interior(field)))
    N₀ = concentrations(model.tracers.N)
    R₀ = concentrations(model.tracers.nitrogen_reservoir)
    inventory₀ = sum(N₀) + sum(R₀)

    @test N₀[1] == 0
    @test R₀[1] ≈ -0.5
    @test N₀[2] == R₀[2] == 0

    for _ in 1:12
        time_step!(model, 0.1)
    end

    N = concentrations(model.tracers.N)
    R = concentrations(model.tracers.nitrogen_reservoir)
    @test all(>=(-1e-12), N)
    @test all(<=(1e-12), R)
    @test R[2] < 0               # debt diffuses into the nutrient-free cell
    @test N[3] < N₀[3]          # some debt reaches and consumes positive N
    @test R[1] > R₀[1]          # debt has left its original cell
    @test sum(N) + sum(R) ≈ inventory₀ atol = 1e-9
end

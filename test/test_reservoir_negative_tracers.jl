using OceanBioME: conserved_tracers
using OceanBioME: reservoir_negative_values
using Oceananigans.TimeSteppers: update_state!
using Adapt: adapt

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
    @test isbitstype(typeof(modifier))
    @test isbitstype(typeof(ScaleNegativeTracers((:N, :P))))
    @test isbitstype(typeof(ClipNegativeTracers(exclude = (:Z,))))
    @test hasproperty(modifier, :tracers) && hasproperty(modifier, :reservoir_names)
    @test hasproperty(ScaleNegativeTracers((:N, :P)), :tracers)
    @test hasproperty(ClipNegativeTracers(exclude = (:Z,)), :exclude)
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
    set!(model.tracers.N, x -> x < 1 ? -0.5 : x < 2 ? 0.0 : 1.0)
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

@testset "Reservoir model compatibility" begin
    grid = RectilinearGrid(architecture, size = (1, 1, 1), extent = (1, 1, 1))
    light = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))

    # LOBSTER can have replicated inorganic carbon, but the organic carbon
    # pools are shared by both conservation groups.
    base = LOBSTER(grid; inorganic_carbon = CarbonateSystem(2), oxygen = Oxygen(),
                   light_attenuation = light)
    modifier = ReservoirNegativeTracers(base.underlying_biogeochemistry;
        repayment_tracers = (nitrogen = :NO₃, carbon1 = :DIC1, carbon2 = :DIC2))
    bgc = Biogeochemistry(base.underlying_biogeochemistry;
                          light_attenuation = light, modifiers = modifier)
    model = NonhydrostaticModel(grid; biogeochemistry = bgc, advection = nothing,
                                auxiliary_fields = (T = ConstantField(15.0), S = ConstantField(35.0)))

    @test modifier.reservoir_names == (:nitrogen_reservoir, :carbon1_reservoir, :carbon2_reservoir)
    @test !(:oxygen_reservoir in keys(model.tracers)) # signed oxygen coefficients are excluded
    @test OceanBioME.modifier_tracers(adapt(CPU(), bgc).modifiers) == modifier.reservoir_names

    set_tracers!(model; NO₃ = 1, NH₄ = 0, P = -0.5, Z = 0.2,
                        DOM = 0.1, sPOM = 0.1, bPOM = 0.1, DIC1 = 0.1, DIC2 = 4,
                        nitrogen_reservoir = 0, carbon1_reservoir = 0, carbon2_reservoir = 0)
    groups = conserved_tracers(bgc)[(:nitrogen, :carbon1, :carbon2)]
    before = group_totals(model, groups)
    update_state!(model)
    after = group_totals(model, groups)

    @test all(name -> isapprox(after[name] + tracer_value(model, Symbol(name, :_reservoir)),
                               before[name]), keys(groups))
    @test all(name -> tracer_value(model, name) >= 0,
              (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM, :DIC1, :DIC2))
    @test tracer_value(model, :carbon1_reservoir) < 0
    @test tracer_value(model, :carbon2_reservoir) ≈ 0

    # PISCES exposes five independent elemental budgets. Verify that inferred
    # coefficients conserve all five even when a shared organic tracer is negative.
    pisces_grid = RectilinearGrid(architecture; topology = (Flat, Flat, Bounded),
                                  size = (10,), extent = (100,))
    warning = "This implementation of PISCES is in early development and has not yet been validated against the operational version"
    pisces = @test_warn warning PISCES(; grid = pisces_grid, light_attenuation = nothing)
    pisces_modifier = ReservoirNegativeTracers(pisces.underlying_biogeochemistry;
        repayment_tracers = (carbon = :DIC, iron = :Fe, phosphate = :PO₄,
                             silicon = :Si, nitrogen = :NO₃))
    @test pisces_modifier.reservoir_names ==
          (:carbon_reservoir, :iron_reservoir, :phosphate_reservoir,
           :silicon_reservoir, :nitrogen_reservoir)
    @test !(:PChl in pisces_modifier.tracers) # excluded from PISCES conservation groups

    c = map(name -> name == :P ? -0.5 : name in (:DIC, :Fe, :PO₄, :Si, :NO₃) ? 2.0 : 0.0,
            pisces_modifier.tracers)
    r = ntuple(_ -> 0.0, length(pisces_modifier.reservoir_names))
    corrected, remaining = reservoir_negative_values(c, r, pisces_modifier.scalefactors,
                                                    pisces_modifier.repayment_indices)
    @test all(>=(0), corrected)
    @test all(<=(0), remaining)
    for g in eachindex(remaining)
        weights = pisces_modifier.scalefactors[g]
        inventory(values) = sum(w * v for (w, v) in zip(weights, values))
        @test inventory(corrected) + remaining[g] ≈ inventory(c)
    end
end

@testset "Reservoir checkpoint restart (CPU)" begin
    # Oceananigans checkpoints the reservoir as an ordinary prognostic tracer.
    # Exercise CPU checkpoint I/O on either test host (CPU or GPU).
    function restart_model()
        grid = RectilinearGrid(CPU(), size = (1, 1, 1), extent = (1, 1, 1))
        base = NPZD(grid; plankton = Abiotic(),
                    detritus = InstantRemineralisationDetritus(), light_attenuation = nothing)
        modifier = ReservoirNegativeTracers(base.underlying_biogeochemistry;
                                           repayment_tracers = (nitrogen = :N,))
        bgc = Biogeochemistry(base.underlying_biogeochemistry; modifiers = modifier)
        model = NonhydrostaticModel(grid; biogeochemistry = bgc, advection = nothing)
        set_tracers!(model; N = 0.1, nitrogen_reservoir = -0.6)
        return model
    end

    mktempdir() do dir
        original = Simulation(restart_model(); Δt = 0.1, stop_iteration = 2)
        original.output_writers[:checkpoint] = Checkpointer(original.model;
            schedule = IterationInterval(2), prefix = "reservoir", dir)
        run!(original)
        @test isfile(joinpath(dir, "reservoir_iteration2.jld2"))

        resumed = Simulation(restart_model(); Δt = 0.1, stop_iteration = 3)
        resumed.output_writers[:checkpoint] = Checkpointer(resumed.model;
            schedule = IterationInterval(2), prefix = "reservoir", dir)
        run!(resumed; pickup = :highest_iteration)

        uninterrupted = Simulation(restart_model(); Δt = 0.1, stop_iteration = 3)
        run!(uninterrupted)
        @test resumed.model.clock.iteration == uninterrupted.model.clock.iteration == 3
        for name in (:N, :nitrogen_reservoir)
            @test tracer_value(resumed.model, name) ≈ tracer_value(uninterrupted.model, name)
        end
        @test tracer_value(resumed.model, :nitrogen_reservoir) < 0
    end
end

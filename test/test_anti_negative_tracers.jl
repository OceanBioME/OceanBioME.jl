using OceanBioME: conserved_tracers
using OceanBioME: anti_negative_values
using Oceananigans.TimeSteppers: update_state!
using Adapt: adapt

@testset "AntiNegativeTracers" begin
    # Pure-cell calculation: a negative total survives as a negative anti-tracer value,
    # and the correction is idempotent when there is nothing left to annihilate it.
    for FT in (Float32, Float64)
        weights = ((one(FT), one(FT)),)
        corrected, anti_tracers = anti_negative_values((FT(0.25), FT(-1)),
                                                      (zero(FT),), weights, (1,))
        @test corrected == (zero(FT), zero(FT))
        @test anti_tracers[1] ≈ FT(-0.75)
        @test sum((FT(0.25), FT(-1))) ≈ sum(corrected) + anti_tracers[1]
        @test anti_negative_values(corrected, anti_tracers, weights, (1,)) == (corrected, anti_tracers)

        # Anti-tracer values use the selected tracer's units, even for a non-unit weight.
        weights = ((FT(2), FT(3)),)
        c, a = (FT(0.5), FT(-0.5)), (FT(-0.1),)
        corrected, anti_tracers = anti_negative_values(c, a, weights, (1,))
        @test corrected == (zero(FT), zero(FT))
        @test anti_tracers[1] ≈ FT(-0.35)
        @test isapprox(sum(weights[1] .* c) + weights[1][1] * a[1],
                       sum(weights[1] .* corrected) + weights[1][1] * anti_tracers[1])
        @test anti_negative_values(corrected, anti_tracers, weights, (1,)) == (corrected, anti_tracers)

        corrected, anti_tracers = anti_negative_values((one(FT), FT(-0.5)), a, weights, (1,))
        @test corrected[1] ≈ FT(0.15)
        @test anti_tracers == (zero(FT),)
    end

    grid = RectilinearGrid(architecture, size = (1, 1, 1), extent = (1, 1, 1))
    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))
    inorganic_carbon = CarbonateSystem()
    base = NPZD(grid; inorganic_carbon, light_attenuation)
    groups = conserved_tracers(base)
    modifier = AntiNegativeTracers(base; annihilation_tracers = (nitrogen = :N, carbon = :DIC))
    @test isbitstype(typeof(modifier))
    @test isbitstype(typeof(ScaleNegativeTracers((:N, :P))))
    @test isbitstype(typeof(ClipNegativeTracers(exclude = (:Z,))))
    @test hasproperty(modifier, :tracers) && hasproperty(modifier, :anti_tracer_names)
    @test hasproperty(ScaleNegativeTracers((:N, :P)), :tracers)
    @test hasproperty(ClipNegativeTracers(exclude = (:Z,)), :exclude)
    biogeochemistry = NPZD(grid; inorganic_carbon, light_attenuation, modifiers = modifier)
    model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing)

    @test (:nitrogen_anti_tracer in keys(model.tracers)) && (:carbon_anti_tracer in keys(model.tracers))
    @test OceanBioME.is_anti_tracer_tendency(modifier, 1, 1, 1, grid, Val(:nitrogen_anti_tracer))
    @test !OceanBioME.is_anti_tracer_tendency(modifier, 1, 1, 1, grid, Val(:N))

    for values in ((N = 2.0, P = -0.5, Z = 1.0, D = 1.0, DIC = 100.0),
                   (N = 0.1, P = -1.0, Z = 0.0, D = 0.0, DIC = 0.0))
        set_tracers!(model; N = 0, P = 0, Z = 0, D = 0, DIC = 0,
                            nitrogen_anti_tracer = 0, carbon_anti_tracer = 0)
        set_tracers!(model; values...)

        before = group_totals(model, groups)
        update_state!(model)
        after = group_totals(model, groups)
        after = (nitrogen = after.nitrogen + tracer_value(model, :nitrogen_anti_tracer),
                carbon = after.carbon + tracer_value(model, :carbon_anti_tracer))

        @test all(map(≈, before, after))
        @test all(name -> tracer_value(model, name) >= 0, (:N, :P, :Z, :D, :DIC))
        @test tracer_value(model, :nitrogen_anti_tracer) <= 0
        @test tracer_value(model, :carbon_anti_tracer) <= 0
    end

    # Annihilation with available nutrient leaves the two conserved inventories unchanged.
    set_tracers!(model; N = 2, P = 0, Z = 0, D = 0, DIC = 10,
                        nitrogen_anti_tracer = -0.5, carbon_anti_tracer = -1)
    update_state!(model)
    @test tracer_value(model, :N) ≈ 1.5
    @test tracer_value(model, :DIC) ≈ 9
    @test tracer_value(model, :nitrogen_anti_tracer) ≈ 0
    @test tracer_value(model, :carbon_anti_tracer) ≈ 0

    @test_throws ArgumentError AntiNegativeTracers(
        (nitrogen = (N = 1, P = 1), carbon = (DIC = 1, P = 6.56));
        annihilation_tracers = (nitrogen = :P, carbon = :DIC))
    @test_throws ArgumentError AntiNegativeTracers(
        (oxygen = (O₂ = 1, P = -7),); annihilation_tracers = (oxygen = :O₂,))
end

@testset "Anti-tracer diffusion and annihilation" begin
    # An abiotic, single-nutrient model isolates Oceananigans transport from
    # biological reactions. The anti-tracer is an ordinary prognostic tracer.
    grid = RectilinearGrid(architecture; size = 16, x = (0, 16),
                           topology = (Bounded, Flat, Flat))
    base = NPZD(grid; plankton = Abiotic(), detritus = InstantRemineralisationDetritus(),
                light_attenuation = nothing)
    modifier = AntiNegativeTracers(base.underlying_biogeochemistry;
                                   annihilation_tracers = (nitrogen = :N,))
    bgc = Biogeochemistry(base.underlying_biogeochemistry; modifiers = modifier)

    closure = HorizontalScalarDiffusivity(ν = 0, κ = (N = 0, nitrogen_anti_tracer = 0.5))
    model = NonhydrostaticModel(grid; biogeochemistry = bgc, closure,
                                momentum_advection = nothing,
                                tracer_advection = (N = Centered(), nitrogen_anti_tracer = nothing))

    # Cell 1 is negative, cell 2 contains no nutrient, and cells 3+ can
    # annihilate the diffusing anti-tracer value. Boundaries are no-flux.
    set!(model.tracers.N, x -> x < 1 ? -0.5 : x < 2 ? 0.0 : 1.0)
    update_state!(model)

    concentrations(field) = vec(Array(interior(field)))
    N₀ = concentrations(model.tracers.N)
    A₀ = concentrations(model.tracers.nitrogen_anti_tracer)
    inventory₀ = sum(N₀) + sum(A₀)

    @test N₀[1] == 0
    @test A₀[1] ≈ -0.5
    @test N₀[2] == A₀[2] == 0

    for _ in 1:12
        time_step!(model, 0.1)
    end

    N = concentrations(model.tracers.N)
    A = concentrations(model.tracers.nitrogen_anti_tracer)
    @test all(>=(-1e-12), N)
    @test all(<=(1e-12), A)
    @test A[2] < 0            # anti-tracer diffuses into the nutrient-free cell
    @test N[3] < N₀[3]       # annihilation reduces N in a neighbouring cell
    @test A[1] > A₀[1]       # anti-tracer value has left its original cell
    @test sum(N) + sum(A) ≈ inventory₀ atol = 1e-9
end

@testset "Anti-tracer model compatibility" begin
    grid = RectilinearGrid(architecture, size = (1, 1, 1), extent = (1, 1, 1))
    light = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))

    # LOBSTER can have replicated inorganic carbon, but the organic carbon
    # pools are shared by both conservation groups.
    base = LOBSTER(grid; inorganic_carbon = CarbonateSystem(2), oxygen = Oxygen(),
                   light_attenuation = light)
    modifier = AntiNegativeTracers(base.underlying_biogeochemistry;
        annihilation_tracers = (nitrogen = :NO₃, carbon1 = :DIC1, carbon2 = :DIC2))
    bgc = Biogeochemistry(base.underlying_biogeochemistry;
                          light_attenuation = light, modifiers = modifier)
    model = NonhydrostaticModel(grid; biogeochemistry = bgc, advection = nothing,
                                auxiliary_fields = (T = ConstantField(15.0), S = ConstantField(35.0)))

    @test modifier.anti_tracer_names == (:nitrogen_anti_tracer, :carbon1_anti_tracer, :carbon2_anti_tracer)
    @test !(:oxygen_anti_tracer in keys(model.tracers)) # groups with negative oxygen weights are excluded
    @test OceanBioME.modifier_tracers(adapt(CPU(), bgc).modifiers) == modifier.anti_tracer_names

    set_tracers!(model; NO₃ = 1, NH₄ = 0, P = -0.5, Z = 0.2,
                        DOM = 0.1, sPOM = 0.1, bPOM = 0.1, DIC1 = 0.1, DIC2 = 4,
                        nitrogen_anti_tracer = 0, carbon1_anti_tracer = 0, carbon2_anti_tracer = 0)
    groups = conserved_tracers(bgc)[(:nitrogen, :carbon1, :carbon2)]
    before = group_totals(model, groups)
    update_state!(model)
    after = group_totals(model, groups)

    @test all(name -> isapprox(after[name] + tracer_value(model, Symbol(name, :_anti_tracer)),
                               before[name]), keys(groups))
    @test all(name -> tracer_value(model, name) >= 0,
              (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM, :DIC1, :DIC2))
    @test tracer_value(model, :carbon1_anti_tracer) < 0
    @test tracer_value(model, :carbon2_anti_tracer) ≈ 0

    # PISCES exposes five independent elemental budgets. Verify that inferred
    # weights conserve all five even when a shared organic tracer is negative.
    pisces_grid = RectilinearGrid(architecture; topology = (Flat, Flat, Bounded),
                                  size = (10,), extent = (100,))
    warning = "This implementation of PISCES is in early development and has not yet been validated against the operational version"
    pisces = @test_warn warning PISCES(; grid = pisces_grid, light_attenuation = nothing)
    pisces_modifier = AntiNegativeTracers(pisces.underlying_biogeochemistry;
        annihilation_tracers = (carbon = :DIC, iron = :Fe, phosphate = :PO₄,
                             silicon = :Si, nitrogen = :NO₃))
    @test pisces_modifier.anti_tracer_names ==
          (:carbon_anti_tracer, :iron_anti_tracer, :phosphate_anti_tracer,
           :silicon_anti_tracer, :nitrogen_anti_tracer)
    @test !(:PChl in pisces_modifier.tracers) # excluded from PISCES conservation groups

    c = map(name -> name == :P ? -0.5 : name in (:DIC, :Fe, :PO₄, :Si, :NO₃) ? 2.0 : 0.0,
            pisces_modifier.tracers)
    a = ntuple(_ -> 0.0, length(pisces_modifier.anti_tracer_names))
    corrected, anti_values = anti_negative_values(c, a, pisces_modifier.scalefactors,
                                                 pisces_modifier.annihilation_indices)
    @test all(>=(0), corrected)
    @test all(<=(0), anti_values)
    for g in eachindex(anti_values)
        weights = pisces_modifier.scalefactors[g]
        inventory(values) = sum(w * v for (w, v) in zip(weights, values))
        alpha = weights[pisces_modifier.annihilation_indices[g]]
        @test inventory(corrected) + alpha * anti_values[g] ≈ inventory(c) + alpha * a[g]
    end
end

@testset "Anti-tracer checkpoint restart (CPU)" begin
    # Oceananigans checkpoints the anti-tracer as an ordinary prognostic tracer.
    # Exercise CPU checkpoint I/O on either test host (CPU or GPU).
    function restart_model()
        grid = RectilinearGrid(CPU(), size = (1, 1, 1), extent = (1, 1, 1))
        base = NPZD(grid; plankton = Abiotic(),
                    detritus = InstantRemineralisationDetritus(), light_attenuation = nothing)
        modifier = AntiNegativeTracers(base.underlying_biogeochemistry;
                                       annihilation_tracers = (nitrogen = :N,))
        bgc = Biogeochemistry(base.underlying_biogeochemistry; modifiers = modifier)
        model = NonhydrostaticModel(grid; biogeochemistry = bgc, advection = nothing)
        set_tracers!(model; N = 0.1, nitrogen_anti_tracer = -0.6)
        return model
    end

    mktempdir() do dir
        original = Simulation(restart_model(); Δt = 0.1, stop_iteration = 2)
        original.output_writers[:checkpoint] = Checkpointer(original.model;
            schedule = IterationInterval(2), prefix = "anti_tracer", dir)
        run!(original)
        @test isfile(joinpath(dir, "anti_tracer_iteration2.jld2"))

        resumed = Simulation(restart_model(); Δt = 0.1, stop_iteration = 3)
        resumed.output_writers[:checkpoint] = Checkpointer(resumed.model;
            schedule = IterationInterval(2), prefix = "anti_tracer", dir)
        run!(resumed; pickup = :highest_iteration)

        uninterrupted = Simulation(restart_model(); Δt = 0.1, stop_iteration = 3)
        run!(uninterrupted)
        @test resumed.model.clock.iteration == uninterrupted.model.clock.iteration == 3
        for name in (:N, :nitrogen_anti_tracer)
            @test tracer_value(resumed.model, name) ≈ tracer_value(uninterrupted.model, name)
        end
        @test tracer_value(resumed.model, :nitrogen_anti_tracer) < 0
    end
end

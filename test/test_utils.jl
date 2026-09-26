include("dependencies_for_runtests.jl")

using OceanBioME: setup_velocity_fields, valid_sinking_velocity_locations, conserved_tracers

using Oceananigans.Fields: AbstractField, CenterField, ConstantField, FunctionField, ZFaceField, location
using Oceananigans.TimeSteppers: update_state!

function test_negative_scaling(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(grid; biogeochemistry = NPZD(grid; scale_negatives = true))

    set!(model, N = 2, P = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)

    run!(simulation)

    N = Array(interior(model.tracers.N))[1, 1, 1]
    P = Array(interior(model.tracers.P))[1, 1, 1]

    return (N ≈ 1) && (P ≈ 0.0)
end

function test_negative_zeroing(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    model = NonhydrostaticModel(grid; biogeochemistry = NPZD(grid; modifiers = ZeroNegativeTracers(; exclude = (:Z, ))))

    set!(model, N = 2, P = -1, Z = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)

    run!(simulation)

    N = Array(interior(model.tracers.N))[1, 1, 1]
    P = Array(interior(model.tracers.P))[1, 1, 1]
    Z = Array(interior(model.tracers.Z))[1, 1, 1]

    return (N ≈ 2) && (P ≈ 0.0) && (Z ≈ -1)
end

@testset "Test negative tracer handeling" begin
    @test test_negative_scaling(architecture)
    @test test_negative_zeroing(architecture)
end

# set the tracers without calling `update_state!`, which is where the negative values are scaled
set_tracers!(model; values...) = [set!(model.tracers[name], value) for (name, value) in pairs(values)]

tracer_value(model, name) = CUDA.@allowscalar model.tracers[name][1, 1, 1]

group_totals(model, groups) = map(group -> sum(scalefactor * tracer_value(model, name) for (name, scalefactor) in pairs(group)), groups)

@testset "ScaleNegativeTracers" begin
    grid = RectilinearGrid(architecture, size = (1, 1, 1), extent = (1, 1, 1))
    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))

    @testset "Groups which share tracers are all conserved" begin
        biogeochemistry = NPZD(grid; inorganic_carbon = CarbonateSystem(), scale_negatives = true, light_attenuation)
        model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing)
        groups = conserved_tracers(biogeochemistry)

        # a negative tracer in both groups, only in the nitrogen group, and only in the carbon group
        for values in ((N = 2, P = -0.5, Z = 1, D = 1, DIC = 2000),
                       (N = -0.5, P = 2, Z = 1, D = 1, DIC = 2000),
                       (N = 2, P = 1, Z = 1, D = 1, DIC = -20))
            set_tracers!(model; values...)
            before = group_totals(model, groups)
            update_state!(model)
            after = group_totals(model, groups)

            @test after.nitrogen ≈ before.nitrogen
            @test after.carbon ≈ before.carbon
            @test all(name -> tracer_value(model, name) ≥ 0, keys(values))
        end

        # the same groups given by hand, with a negative value that makes the nitrogen group scale the plankton
        modifiers = ScaleNegativeTracers((nitrogen = (:N, :P, :Z, :D), carbon = (DIC = 1, P = 6.56, Z = 6.56, D = 6.56)))
        biogeochemistry = NPZD(grid; inorganic_carbon = CarbonateSystem(), modifiers, light_attenuation)
        model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing)

        set_tracers!(model; N = -0.5, P = 2, Z = 1, D = 1, DIC = 2000)
        update_state!(model)

        @test tracer_value(model, :N) == 0
        @test tracer_value(model, :P) + tracer_value(model, :Z) + tracer_value(model, :D) ≈ 3.5
        @test tracer_value(model, :DIC) + 6.56 * 3.5 ≈ 2000 + 6.56 * 4
    end

    @testset "Scale factors" begin
        modifiers = ScaleNegativeTracers((:P, :Z, :N); scalefactors = (1, 1, 2))
        model = NonhydrostaticModel(grid; biogeochemistry = NPZD(grid; modifiers, light_attenuation), advection = nothing)

        set_tracers!(model; N = 1, P = -1, Z = 1)
        update_state!(model)

        @test tracer_value(model, :P) == 0
        @test tracer_value(model, :Z) + 2 * tracer_value(model, :N) ≈ 2

        @test_throws ArgumentError ScaleNegativeTracers((O₂ = 1, P = -7))
    end

    @testset "Replicated inorganic carbon" begin
        for inorganic_carbon in (CarbonateSystem(2), ExplicitCalciumCarbonate(grid; replicates = 2))
            biogeochemistry = LOBSTER(grid; inorganic_carbon, scale_negatives = true, light_attenuation)
            model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing,
                                        auxiliary_fields = (T = ConstantField(15.0), S = ConstantField(35.0)))
            groups = conserved_tracers(biogeochemistry)

            @test keys(groups) == (:nitrogen, :carbon1, :carbon2)

            for (DIC1, DIC2) in ((2000, 2000), (2000, 2100))
                set_tracers!(model; NO₃ = 10, NH₄ = 0.1, P = -0.5, Z = 0.5, DOM = 0.2, sPOM = 0.1, bPOM = 0.1, DIC1, DIC2)
                inorganic_carbon isa ExplicitCalciumCarbonate && set_tracers!(model; CaCO₃1 = 1, CaCO₃2 = 1)

                before = group_totals(model, groups)
                update_state!(model)
                after = group_totals(model, groups)

                @test all(map(≈, after, before))

                # identical replicates stay identical
                @test (tracer_value(model, :DIC1) == tracer_value(model, :DIC2)) == (DIC1 == DIC2)
            end
        end
    end

    @testset "Oxygen is not scaled" begin
        biogeochemistry = LOBSTER(grid; oxygen = Oxygen(), scale_negatives = true, light_attenuation)
        model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing)

        # anoxic water with more organic matter than its nitrate can oxidise, but nothing negative
        set_tracers!(model; NO₃ = 1, NH₄ = 0, P = 0, Z = 0, DOM = 1, sPOM = 0, bPOM = 0, O₂ = 0)
        update_state!(model)

        @test tracer_value(model, :NO₃) == 1
        @test tracer_value(model, :DOM) == 1

        # negative oxygen does not take nitrogen with it
        set_tracers!(model; NO₃ = 30, NH₄ = 0.1, DOM = 2, sPOM = 1, bPOM = 0.5, O₂ = -1)
        update_state!(model)

        @test sum(name -> tracer_value(model, name), (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM)) ≈ 33.6
        @test tracer_value(model, :O₂) == -1
    end

    @testset "Halos are filled after scaling" begin
        grid = RectilinearGrid(architecture, size = 4, z = (-4, 0), topology = (Flat, Flat, Bounded))
        model = NonhydrostaticModel(grid; biogeochemistry = NPZD(grid; scale_negatives = true),
                                    advection = nothing, closure = ScalarDiffusivity(κ = 0.25),
                                    timestepper = :QuasiAdamsBashforth2)

        set!(model.tracers.N, 1)
        set!(model.tracers.P, z -> z > -3 ? 1 : -0.5)
        update_state!(model)

        # the halo below a no-flux bottom holds a copy of the bottom cell
        @test CUDA.@allowscalar model.tracers.P[1, 1, 0] == model.tracers.P[1, 1, 1] == 0

        # so diffusion between the no-flux boundaries does not change the total
        set!(model.tracers.N, 1)
        set!(model.tracers.P, z -> z > -3 ? 1 : -0.5)
        time_step!(model, 1)

        @test isapprox(sum(Array(interior(model.tracers.P))), 3, atol = 1e-4)
    end
end

scalar_sinking_speeds = (A = 1, B = 1.0)

grid = RectilinearGrid(architecture, size = (1, 1, 10), extent = (1, 1, 10))

field_sinking_speeds = (C = ConstantField(1), D = FunctionField{Center, Center, Face}((x, y, z) -> z, grid), E = ZFaceField(grid))

sinking_speeds = merge(scalar_sinking_speeds, field_sinking_speeds)

@testset "Test sinking velocity setup" begin
    sinking_velocities = @test_nowarn setup_velocity_fields(sinking_speeds, grid, true)

    @test all(map(w -> isa(w, AbstractField) & (location(w) in valid_sinking_velocity_locations), values(sinking_velocities)))

    sinking_velocities =
        @test_warn ("The sinking velocity provided for C is a field and therefore `open_bottom=false` can't be enforced automatically",
                    "The sinking velocity provided for D is a field and therefore `open_bottom=false` can't be enforced automatically",
                    "The sinking velocity provided for E is a field and therefore `open_bottom=false` can't be enforced automatically") setup_velocity_fields(sinking_speeds, grid, false)

    @test all([isa(w, AbstractField) & (location(w) in valid_sinking_velocity_locations) for w in values(sinking_velocities)])

    @test all(map(w -> Array(interior(w, 1, 1, 1)) .== 0, sinking_velocities[(:A, :B)]))

    @test_warn "The location of the sinking velocity field provided for X is incorrect, it should be (Center, Center, Face)" setup_velocity_fields((X = CenterField(grid), ), grid, true)
end

using Oceananigans.Grids: znode, Center
using OceanBioME.Models.NutrientsPlanktonDetritusModels: InstantRemineralisationDetritus

@testset "IronDustDeposition" begin
    grid = RectilinearGrid(architecture; size=(1, 1, 2), extent=(1, 1, 4))
    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))

    dust_flux = 1e-8

    biogeochemistry = NutrientsPlanktonDetritus(grid;
                                                plankton = Abiotic(),
                                                nutrients = Nutrients(; phosphate = OceanBioME.PO₄,
                                                                        iron = SimpleIron(scavenging_rate = 0.0)),
                                                detritus = InstantRemineralisationDetritus(),
                                                light_attenuation)

    iron_forcing = IronDustDepositionForcing(dust_flux)

    model = NonhydrostaticModel(grid;
                                advection = nothing,
                                biogeochemistry,
                                forcing = (; Fe = iron_forcing))

    set!(model, Fe = 0.0, PO₄ = 1.0)
    time_step!(model, 1.0)
    @test CUDA.@allowscalar model.tracers.Fe[1, 1, 2] > 0
    @test CUDA.@allowscalar model.tracers.Fe[1, 1, 2] > model.tracers.Fe[1, 1, 1]

    f_iron = 0.035
    M_Fe = 55.845e-6
    λ = 400.0

    # two-component MARBL formulation (default)
    γ = 0.98
    λₕ = 1.2e6
    dep = IronDustDeposition(dust_flux)
    z = znode(1, 1, 2, grid, Center(), Center(), Center())
    expected = dust_flux * f_iron / M_Fe * ((1 - γ) / λ * exp(z / λ) + γ / λₕ * exp(z / λₕ))
    @test dep(1, 1, 2, grid, nothing, nothing) ≈ expected

    # single exponential
    dep_single = IronDustDeposition(dust_flux; hard_fraction = 0.0)
    expected_single = dust_flux * f_iron / M_Fe / λ * exp(z / λ)
    @test dep_single(1, 1, 2, grid, nothing, nothing) ≈ expected_single
    @test dep_single(1, 1, 2, grid, nothing, nothing) > dep(1, 1, 2, grid, nothing, nothing)
end

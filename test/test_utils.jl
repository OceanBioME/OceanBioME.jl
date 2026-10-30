include("dependencies_for_runtests.jl")

using OceanBioME: setup_velocity_fields, valid_sinking_velocity_locations, conserved_tracers

using Oceananigans.Fields: AbstractField, CenterField, ConstantField, FunctionField, ZFaceField, location
using Oceananigans.TimeSteppers: update_state!

scalar_value(field) = CUDA.@allowscalar field[1, 1, 1]

struct ContinuousNegativeTracerTestBGC <: Oceananigans.Biogeochemistry.AbstractContinuousFormBiogeochemistry end

Oceananigans.Biogeochemistry.required_biogeochemical_tracers(::ContinuousNegativeTracerTestBGC) = (:A, :T)

@inline (::ContinuousNegativeTracerTestBGC)(::Val{:A}, x, y, z, t, A, T, auxiliary) = A + T + auxiliary

function test_negative_scaling(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    bgc = NPZD(grid; negative_tracers = ScaleNegativeTracers())
    model = NonhydrostaticModel(grid; biogeochemistry = bgc)

    set!(model, N = 2, P = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)
    run!(simulation)

    N = scalar_value(model.tracers.N)
    P = scalar_value(model.tracers.P)

    return (N ≈ 1) && (P ≈ 0.0)
end

function test_negative_clipping(arch)
    grid = RectilinearGrid(arch, size = (1, 1, 1), extent = (1, 1, 1))

    bgc = NPZD(grid; negative_tracers = ClipNegativeTracers(; exclude = (:Z, )))
    model = NonhydrostaticModel(grid; biogeochemistry = bgc)

    set!(model, N = 2, P = -1, Z = -1)

    simulation = Simulation(model, Δt = 1e-10, stop_iteration = 1)
    run!(simulation)

    N = scalar_value(model.tracers.N)
    P = scalar_value(model.tracers.P)
    Z = scalar_value(model.tracers.Z)

    return (N ≈ 2) && (P ≈ 0.0) && (Z ≈ -1)
end

@testset "Test negative tracer handling" begin
    @test test_negative_scaling(architecture)
    @test test_negative_clipping(architecture)

    grid = RectilinearGrid(architecture, size = (1, 1, 1), extent = (1, 1, 1))
    @test NPZD(grid; scale_negatives = true).negative_tracers isa ScaleNegativeTracers
    @test_throws ArgumentError NPZD(grid; negative_tracers = IgnoreNegativeTracerValues(), scale_negatives = true)

    raw = Biogeochemistry(ContinuousNegativeTracerTestBGC())
    safe = Biogeochemistry(ContinuousNegativeTracerTestBGC(); negative_tracers = IgnoreNegativeTracerValues())

    @test raw(Val(:A), 0, 0, 0, 0, -1.0, -2.0, -3.0) == -6.0
    @test safe(Val(:A), 0, 0, 0, 0, -1.0, -2.0, -3.0) == -5.0
    composed = Biogeochemistry(ContinuousNegativeTracerTestBGC();
                                negative_tracers = (ClipNegativeTracers(), (IgnoreNegativeTracerValues(),)))
    @test composed(Val(:A), 0, 0, 0, 0, -1.0, -2.0, -3.0) == -5.0

    light = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(10.0))
    reference = NPZD(grid; light_attenuation = light)
    safe_npzd = NPZD(grid; light_attenuation = light, negative_tracers = IgnoreNegativeTracerValues())
    reference_model = NonhydrostaticModel(grid; biogeochemistry = reference)
    safe_model = NonhydrostaticModel(grid; biogeochemistry = safe_npzd)

    common = (; N = 2.0, Z = 0.1, D = 0.1, T = -1.0)
    set!(reference_model; common..., P = 0.0)
    set!(safe_model; common..., P = -1e-3)

    tendency_reference = CUDA.@allowscalar reference(1, 1, 1, grid, Val(:P), reference_model.clock, reference_model.tracers)
    tendency_safe = CUDA.@allowscalar safe_npzd(1, 1, 1, grid, Val(:P), safe_model.clock, safe_model.tracers)

    @test tendency_safe ≈ tendency_reference

    simulation = Simulation(safe_model, Δt = 1e-10, stop_iteration = 1)
    run!(simulation)
    @test scalar_value(safe_model.tracers.P) ≈ -1e-3

    for light_model in (TwoBandPhotosyntheticallyActiveRadiation,
                        MultiBandPhotosyntheticallyActiveRadiation)
        reference = NPZD(grid; light_attenuation = light_model(grid, 100.0))
        safe_npzd = NPZD(grid;
                         light_attenuation = light_model(grid, 100.0),
                         negative_tracers = IgnoreNegativeTracerValues())
        reference_model = NonhydrostaticModel(grid; biogeochemistry = reference)
        safe_model = NonhydrostaticModel(grid; biogeochemistry = safe_npzd)

        set!(reference_model; common..., P = 0.0)
        set!(safe_model; common..., P = -1e-3)

        Oceananigans.Biogeochemistry.update_biogeochemical_state!(reference_model, reference.light_attenuation)
        Oceananigans.Biogeochemistry.update_biogeochemical_state!(safe_model, safe_npzd.light_attenuation)

        reference_light = Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(reference.light_attenuation)
        safe_light = Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(safe_npzd.light_attenuation)

        @test all(name -> scalar_value(safe_light[name]) ≈ scalar_value(reference_light[name]), keys(reference_light))
    end
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

        # A negative conserved total cannot be made nonnegative while conserving it.
        set_tracers!(model; N = 0.1, P = -1, Z = 0)
        update_state!(model)
        @test isnan(tracer_value(model, :N))

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

    @testset "Temperature is not scaled" begin
        biogeochemistry = NPZD(grid; scale_negatives = true, light_attenuation)
        model = NonhydrostaticModel(grid; biogeochemistry, advection = nothing)

        @test keys(conserved_tracers(biogeochemistry).nitrogen) == (:N, :P, :Z, :D)

        set_tracers!(model; N = 5, P = 0.1, Z = 0.5, D = 0.5, T = -1)
        update_state!(model)

        @test tracer_value(model, :T) == -1
        @test tracer_value(model, :N) == 5
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

@testset "A closed bottom is closed at the seafloor of an immersed grid" begin
    # Column 1 has `Nz_wet` wet cells above the seafloor and column 2 is land. With a closed bottom the sinking velocity
    # above the seafloor must be the same as on a grid whose bottom is the seafloor, and zero at and below the seafloor
    Nz, Nz_wet = 10, 4
    kᵇ = Nz - Nz_wet + 1 # the lowest wet cell, whose bottom face is the seafloor

    seafloor(x, y) = ifelse(x < 1, -Nz_wet, 0)

    for FT in (Float64, Float32), bottom in (GridFittedBottom(seafloor), PartialCellBottom(seafloor))
        grid = ImmersedBoundaryGrid(RectilinearGrid(architecture, FT; size = (2, 1, Nz), extent = (2, 1, Nz)), bottom)
        shallow_grid = RectilinearGrid(architecture, FT; size = (2, 1, Nz_wet), extent = (2, 1, Nz_wet))

        w_closed = Array(interior(setup_velocity_fields((; A = 1), grid, false).A))
        w_shallow = Array(interior(setup_velocity_fields((; A = 1), shallow_grid, false).A))
        w_open = Array(interior(setup_velocity_fields((; A = 1), grid, true).A))

        @test w_closed[1, 1, kᵇ:end] == w_shallow[1, 1, :]
        @test all(w_closed[1, 1, 1:kᵇ] .== 0)
        @test all(w_closed[2, 1, :] .== 0)
        @test all(w_open[1, 1, kᵇ:Nz] .== -1) # an open bottom still lets tracers sink through the seafloor
    end
end

@testset "Sinking tracers are conserved with a closed bottom on an immersed grid" begin
    # Detritus sinks and remineralises into N, so with a closed bottom N + D can't change, and the
    # water above the seafloor must evolve like it does on a grid whose bottom is the seafloor
    deep_grid = RectilinearGrid(architecture; size = 8, z = (-8, 0), topology = (Flat, Flat, Bounded))
    immersed_grid = ImmersedBoundaryGrid(deep_grid, GridFittedBottom(-4))
    shallow_grid = RectilinearGrid(architecture; size = 4, z = (-4, 0), topology = (Flat, Flat, Bounded))

    immersed, shallow = map((immersed_grid, shallow_grid)) do grid
        detritus = Detritus(grid; sinking_speed = 1, remineralisation_rate = 0.1, open_bottom = false)
        biogeochemistry = NutrientsPlanktonDetritus(grid; nutrients = Nutrients(; nitrogen = OceanBioME.N), detritus)

        # upwinding keeps the detritus piling up on the seafloor positive
        model = HydrostaticFreeSurfaceModel(grid; biogeochemistry, buoyancy = nothing, tracers = (),
                                            tracer_advection = UpwindBiased(order = 1))
        set!(model, D = 1)

        for _ in 1:100
            time_step!(model, 0.1)
        end

        return model
    end

    total_nitrogen = sum(CUDA.@allowscalar first(Field(Integral(immersed.tracers[c]))) for c in (:N, :D))

    @test total_nitrogen ≈ 4 # D = 1 in the 4 m of water above the seafloor
    @test Array(interior(immersed.tracers.D, 1, 1, 5:8)) ≈ Array(interior(shallow.tracers.D, 1, 1, :))
    @test Array(interior(immersed.tracers.N, 1, 1, 5:8)) ≈ Array(interior(shallow.tracers.N, 1, 1, :))
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

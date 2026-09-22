include("dependencies_for_runtests.jl")

using OceanBioME: setup_velocity_fields, valid_sinking_velocity_locations

using Oceananigans.Fields: AbstractField, CenterField, ConstantField, FunctionField, ZFaceField, location

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

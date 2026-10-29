include("dependencies_for_runtests.jl")

using OceanBioME: SimpleCaCO3Precipitation
using Oceananigans

function sum_tracer(model, name)
    return sum(Array(interior(model.tracers[name])))
end

# ---------------------------------------------------------------------------
# BoxModelGrid regression: calcite saturation update should not error
# ---------------------------------------------------------------------------

@testset "SimpleCaCO3Precipitation – BoxModelGrid pressure fallback" begin
    grid = BoxModelGrid()
    bgc  = SimpleCaCO3Precipitation(; grid, sinking_speed = 0)
    model = BoxModel(; biogeochemistry = bgc, clock = Clock(time = 0.0))

    set!(model; DIC = 2000.0, Alk = 2300.0, CaCO₃ = 0.0, T = 25.0, S = 35.0)

    # Should update Ω and take a timestep without trying `abs(::Nothing)`.
    @test begin
        time_step!(model, 60.0)
        true
    end
end

# ---------------------------------------------------------------------------
# Basic instantiation tests
# ---------------------------------------------------------------------------

@testset "SimpleCaCO3Precipitation – instantiation" begin
    grid = RectilinearGrid(architecture; size=(3, 3, 6), extent=(100, 100, 200),
                           topology=(Periodic, Periodic, Bounded))
    bgc = SimpleCaCO3Precipitation(; grid)

    required = Oceananigans.Biogeochemistry.required_biogeochemical_tracers(bgc)
    @test :DIC             ∈ required
    @test :Alk             ∈ required
    @test Symbol("CaCO₃") ∈ required
    @test :T               ∈ required
    @test :S               ∈ required

    aux = Oceananigans.Biogeochemistry.required_biogeochemical_auxiliary_fields(bgc)
    @test :Ω ∈ aux

    ubgc = bgc.underlying_biogeochemistry
    @test ubgc.precipitation_rate_constant isa Float64
    @test ubgc.precipitation_order         isa Float64
    @test ubgc.dissolution_rate_constant   isa Float64
    @test ubgc.dissolution_order           isa Float64
end

# ---------------------------------------------------------------------------
# Float32 parameter types propagate correctly
# ---------------------------------------------------------------------------

@testset "SimpleCaCO3Precipitation – Float32" begin
    grid = RectilinearGrid(architecture, Float32; size=(2, 2, 4), extent=(10, 10, 40),
                           topology=(Periodic, Periodic, Bounded))
    bgc  = SimpleCaCO3Precipitation(; grid)
    ubgc = bgc.underlying_biogeochemistry

    @test ubgc.precipitation_rate_constant isa Float32
    @test ubgc.precipitation_order         isa Float32
    @test ubgc.dissolution_rate_constant   isa Float32
    @test ubgc.dissolution_order           isa Float32
end

# ---------------------------------------------------------------------------
# Sinking velocity registered for CaCO₃, zero for dissolved tracers
# ---------------------------------------------------------------------------

@testset "SimpleCaCO3Precipitation – sinking" begin
    grid = RectilinearGrid(architecture; size=(2, 2, 4), extent=(10, 10, 40),
                           topology=(Periodic, Periodic, Bounded))
    bgc  = SimpleCaCO3Precipitation(; grid)
    ubgc = bgc.underlying_biogeochemistry

    drift = Oceananigans.Biogeochemistry.biogeochemical_drift_velocity(ubgc, Val(Symbol("CaCO₃")))
    @test haskey(drift, :w)

    drift_DIC = Oceananigans.Biogeochemistry.biogeochemical_drift_velocity(ubgc, Val(:DIC))
    @test drift_DIC.w isa Oceananigans.Fields.ZeroField
end

# ---------------------------------------------------------------------------
# Zero initial conditions → zero tendencies
# ---------------------------------------------------------------------------

@testset "SimpleCaCO3Precipitation – zero tracers give zero tendencies" begin
    grid  = RectilinearGrid(architecture; size=(2, 2, 4), extent=(10, 10, 40),
                            topology=(Periodic, Periodic, Bounded))
    model = NonhydrostaticModel(grid; biogeochemistry = SimpleCaCO3Precipitation(; grid))

    time_step!(model, 1.0)

    @test all(Array(interior(model.tracers.DIC))             .== 0)
    @test all(Array(interior(model.tracers.Alk))             .== 0)
    @test all(Array(interior(model.tracers[Symbol("CaCO₃")])).== 0)
end

# ---------------------------------------------------------------------------
# Supersaturation drives precipitation (DIC ↓, Alk ↓, CaCO₃ ↑)
# and carbon budget is conserved when there is no sinking
# ---------------------------------------------------------------------------

@testset "SimpleCaCO3Precipitation – supersaturation precipitation" begin
    grid  = RectilinearGrid(architecture; size=(2, 2, 2), extent=(10, 10, 20),
                            topology=(Periodic, Periodic, Bounded))
    # Disable sinking to cleanly check carbon budget
    bgc   = SimpleCaCO3Precipitation(; grid, sinking_speed = 0)
    model = NonhydrostaticModel(grid; biogeochemistry = bgc)

    # Typical surface seawater: Ω ≈ 5 (supersaturated)
    model.tracers.DIC              .= 2000.0
    model.tracers.Alk              .= 2300.0
    model.tracers[Symbol("CaCO₃")] .= 0.0
    model.tracers.T                .= 25.0
    model.tracers.S                .= 35.0

    DIC0   = sum_tracer(model, :DIC)
    Alk0   = sum_tracer(model, :Alk)
    CaCO30 = sum_tracer(model, Symbol("CaCO₃"))

    for _ in 1:10
        time_step!(model, 3600.0)
    end

    DIC1   = sum_tracer(model, :DIC)
    Alk1   = sum_tracer(model, :Alk)
    CaCO31 = sum_tracer(model, Symbol("CaCO₃"))

    # Precipitation: DIC ↓, Alk ↓, CaCO₃ ↑
    @test DIC1   < DIC0
    @test Alk1   < Alk0
    @test CaCO31 > CaCO30

    # Stoichiometry: ΔAlk = 2 × ΔDIC (both negative)
    ΔDIC = DIC1 - DIC0
    ΔAlk = Alk1 - Alk0
    @test isapprox(ΔAlk / ΔDIC, 2.0; rtol = 1e-6)

    # Carbon budget: DIC + CaCO₃ conserved (no sinking, no external fluxes)
    ΔC_total = (DIC1 + CaCO31) - (DIC0 + CaCO30)
    @test isapprox(ΔC_total, 0.0; atol = 1e-8)
end

# ---------------------------------------------------------------------------
# Undersaturation drives dissolution (DIC ↑, Alk ↑, CaCO₃ ↓)
# ---------------------------------------------------------------------------

@testset "SimpleCaCO3Precipitation – undersaturation dissolution" begin
    grid  = RectilinearGrid(architecture; size=(2, 2, 2), extent=(10, 10, 20),
                            topology=(Periodic, Periodic, Bounded))
    bgc   = SimpleCaCO3Precipitation(; grid, sinking_speed = 0)
    model = NonhydrostaticModel(grid; biogeochemistry = bgc)

    # Deep ocean conditions: DIC=2200, Alk=2000, T=2, S=35 → Ω ≈ 0.23
    model.tracers.DIC              .= 2200.0
    model.tracers.Alk              .= 2000.0
    model.tracers[Symbol("CaCO₃")] .= 100.0
    model.tracers.T                .= 2.0
    model.tracers.S                .= 35.0

    DIC0   = sum_tracer(model, :DIC)
    Alk0   = sum_tracer(model, :Alk)
    CaCO30 = sum_tracer(model, Symbol("CaCO₃"))

    for _ in 1:10
        time_step!(model, 3600.0)
    end

    DIC1   = sum_tracer(model, :DIC)
    Alk1   = sum_tracer(model, :Alk)
    CaCO31 = sum_tracer(model, Symbol("CaCO₃"))

    # Dissolution: DIC ↑, Alk ↑, CaCO₃ ↓
    @test DIC1   > DIC0
    @test Alk1   > Alk0
    @test CaCO31 < CaCO30

    # Stoichiometry: ΔAlk = 2 × ΔDIC (both positive)
    ΔDIC = DIC1 - DIC0
    ΔAlk = Alk1 - Alk0
    @test isapprox(ΔAlk / ΔDIC, 2.0; rtol = 1e-6)

    # Carbon budget conserved
    ΔC_total = (DIC1 + CaCO31) - (DIC0 + CaCO30)
    @test isapprox(ΔC_total, 0.0; atol = 1e-8)
end

# ---------------------------------------------------------------------------
# CaCO₃ stays non-negative under aggressive dissolution
# ---------------------------------------------------------------------------

@testset "SimpleCaCO3Precipitation – CaCO₃ non-negativity" begin
    grid  = RectilinearGrid(architecture; size=(2, 2, 2), extent=(10, 10, 20),
                            topology=(Periodic, Periodic, Bounded))
    bgc   = SimpleCaCO3Precipitation(; grid,
                                       dissolution_rate_constant = 1e-2 / day,
                                       sinking_speed = 0)
    model = NonhydrostaticModel(grid; biogeochemistry = bgc)

    model.tracers.DIC              .= 2200.0
    model.tracers.Alk              .= 2000.0
    model.tracers[Symbol("CaCO₃")] .= 1.0
    model.tracers.T                .= 2.0
    model.tracers.S                .= 35.0

    for _ in 1:1000
        time_step!(model, 3600.0)
    end

    @test all(Array(interior(model.tracers[Symbol("CaCO₃")])) .>= -1e-10)
end

# Prototype tests; intended to be folded into OceanBioME's existing test suite.

using Test
using OceanBioME
using Oceananigans
using Oceananigans.Biogeochemistry: AbstractContinuousFormBiogeochemistry
import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

struct TwoTracerBGC <: AbstractContinuousFormBiogeochemistry end
required_biogeochemical_tracers(::TwoTracerBGC) = (:A, :B)
@inline (::TwoTracerBGC)(::Val{:A}, x, y, z, t, A, B) = A - B
@inline (::TwoTracerBGC)(::Val{:B}, x, y, z, t, A, B) = A + B

@testset "Nonnegative biological state: continuous form" begin
    raw = Biogeochemistry(TwoTracerBGC())
    safe = Biogeochemistry(TwoTracerBGC(); biological_state = NonnegativeBiologicalState())

    @test raw(Val(:A), 0, 0, 0, 0, -1.0, 2.0) == -3.0
    @test safe(Val(:A), 0, 0, 0, 0, -1.0, 2.0) == -2.0
    @test safe(Val(:B), 0, 0, 0, 0, -1.0, 2.0) == 2.0

    # Positive states are exactly unchanged.
    @test safe(Val(:A), 0, 0, 0, 0, 1.0, 2.0) == raw(Val(:A), 0, 0, 0, 0, 1.0, 2.0)
end

@testset "Nonnegative biological state: discrete LOBSTER" begin
    grid = RectilinearGrid(CPU(); size = (1, 1, 1), extent = (1, 1, 1))

    safe = LOBSTER(grid; biological_state = NonnegativeBiologicalState())
    reference = LOBSTER(grid)

    safe_model = NonhydrostaticModel(; grid, biogeochemistry = safe, advection = nothing)
    ref_model = NonhydrostaticModel(; grid, biogeochemistry = reference, advection = nothing)

    # Same physical state except safe keeps a negative phytoplankton value,
    # while the reference receives the value biology is meant to see.
    common = (; NO₃ = 5.0, NH₄ = 0.2, Z = 0.1, DOM = 0.1, sPOM = 0.1, bPOM = 0.1)
    set!(safe_model; common..., P = -1e-3)
    set!(ref_model; common..., P = 0.0)

    aux_safe = Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(safe)
    aux_ref = Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(reference)

    # Wrapper-level tendency should see P=0 without mutating P itself.
    tendency_safe = safe(1, 1, 1, grid, Val(:P), safe_model.clock, safe_model.tracers)
    tendency_ref = reference(1, 1, 1, grid, Val(:P), ref_model.clock, ref_model.tracers)

    @test tendency_safe ≈ tendency_ref
    @test safe_model.tracers.P[1, 1, 1] == -1e-3
end

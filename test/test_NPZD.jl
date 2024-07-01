include("dependencies_for_runtests.jl")

using OceanBioME: NutrientPhytoplanktonZooplanktonDetritus
using Oceananigans

function test_NPZD(grid, sinking, open_bottom)
    PAR = CenterField(grid)

    if sinking
        model = NonhydrostaticModel(;grid,
                                     biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(;grid, open_bottom),
                                     auxiliary_fields = (; PAR))
    else
        model = NonhydrostaticModel(;grid,
                                     biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(;grid, sinking_speeds = NamedTuple()),
                                     auxiliary_fields = (; PAR))
    end

    # correct tracers and auxiliary fields have been setup, and order has not changed
    required_tracers = (:N, :P, :Z, :D, :T)

    @test Oceananigans.Biogeochemistry.required_biogeochemical_tracers(model.biogeochemistry) == required_tracers
    @test all(tracer ∈ keys(model.tracers) for tracer in required_tracers)
    @test :PAR ∈ keys(model.auxiliary_fields)

    # checks model works with zero values
    time_step!(model, 1.0)

    # and that they all return zero
    @test all([all(Array(interior(values)) .== 0) for values in values(model.tracers)]) 

    # mass conservation
    model.tracers.N .= rand()
    model.tracers.P .= rand()
    model.tracers.Z .= rand()
    model.tracers.D .= rand()

    ΣN₀ = sum(Array(interior(model.tracers.N))) + 
          sum(Array(interior(model.tracers.P))) + 
          sum(Array(interior(model.tracers.Z))) + 
          sum(Array(interior(model.tracers.D)))

    for n in 1:1000
        time_step!(model, 1.0)
    end

    ΣN₁ = sum(Array(interior(model.tracers.N))) + 
          sum(Array(interior(model.tracers.P))) + 
          sum(Array(interior(model.tracers.Z))) + 
          sum(Array(interior(model.tracers.D)))

    @test ΣN₀ ≈ ΣN₁ # guess this should actually fail with a high enough accuracy when sinking is on with an open bottom

    return nothing
end

grid = RectilinearGrid(architecture; size=(3, 3, 6), extent=(1, 1, 2))

for sinking = (false, true), open_bottom = (false, true)
    if !(sinking && open_bottom) # no sinking is the same with and without open bottom
        @info "Testing on $(typeof(architecture)) with sinking $(sinking ? :✅ : :❌), open bottom $(open_bottom ? :✅ : :❌))"
        @testset "NPZD $architecture, $sinking, $open_bottom" begin
            test_NPZD(grid, sinking, open_bottom)
        end
    end
end

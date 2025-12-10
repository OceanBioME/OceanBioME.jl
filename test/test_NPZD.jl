include("dependencies_for_runtests.jl")

using OceanBioME: NutrientPhytoplanktonZooplanktonDetritus
using Oceananigans

function test_NPZD(grid, sinking, open_bottom)
    
    if sinking
        model = NonhydrostaticModel(;grid,
                                     biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(;grid, open_bottom))
    else
        model = NonhydrostaticModel(;grid,
                                     biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(;grid, sinking_speeds = NamedTuple()))
    end

    # correct tracers and auxiliary fields have been setup, and order has not changed
    required_tracers = (:N, :P, :Z, :D, :T)

    @test Oceananigans.Biogeochemistry.required_biogeochemical_tracers(model.biogeochemistry) == required_tracers
    @test all(tracer ∈ keys(model.tracers) for tracer in required_tracers)
    @test :PAR ∈ keys(Oceananigans.Biogeochemistry.biogeochemical_auxiliary_fields(model.biogeochemistry))

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

@testset "Float32 NPZD" begin
    grid = RectilinearGrid(architecture, Float32; size=(3, 3, 10), extent=(10, 10, 200))
    bgc = NutrientPhytoplanktonZooplanktonDetritus(; grid)

    ubgc = bgc.underlying_biogeochemistry
    @test ubgc.initial_photosynthetic_slope isa Float32
    @test ubgc.base_maximum_growth isa Float32
    @test ubgc.nutrient_half_saturation isa Float32
    @test ubgc.base_respiration_rate isa Float32
    @test ubgc.phyto_base_mortality_rate isa Float32
    @test ubgc.maximum_grazing_rate isa Float32
    @test ubgc.grazing_half_saturation isa Float32
    @test ubgc.assimulation_efficiency isa Float32
    @test ubgc.base_excretion_rate isa Float32
    @test ubgc.zoo_base_mortality_rate isa Float32
    @test ubgc.remineralization_rate isa Float32

    par = bgc.light_attenuation
    @test par.water_red_attenuation isa Float32
    @test par.water_blue_attenuation isa Float32
    @test par.chlorophyll_red_attenuation isa Float32
    @test par.chlorophyll_blue_attenuation isa Float32
    @test par.chlorophyll_red_exponent isa Float32
    @test par.chlorophyll_blue_exponent isa Float32
    @test par.pigment_ratio isa Float32
    @test par.phytoplankton_chlorophyll_ratio isa Float32
end

using Test
using OceanBioME: LOBSTER
using OceanBioME.Light: TwoBandPhotosyntheticallyActiveRatiationField
using Oceananigans

function test_LOBSTER(grid, carbonates, oxygen, sinking, open_bottom)
    if sinking
        model = NonhydrostaticModel(;grid,
                                    biogeochemistry = LOBSTER(;grid, carbonates, oxygen, open_bottom))
    else
        model = NonhydrostaticModel(;grid,
                                    biogeochemistry = LOBSTER(;grid, carbonates, oxygen, sinking_velocities = NamedTuple()))
    end

    # correct tracers and auxiliary fields have been setup, and order has not changed
    required_tracers = @show !(carbonates && oxygen) ? (:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM) : 
                            (carbonates && !oxygen ? (:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM, :DIC, :ALK) :
                            (oxygen && !carbonates ? (:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM, :OXY) :
                            (:NO₃, :NH₄, :P, :Z, :D, :DD, :Dᶜ, :DDᶜ, :DOM, :DIC, :ALK, :OXY)))

    @test Oceananigans.Biogeochemistry.required_biogeochemical_tracers(model.biogeochemistry) == required_tracers
    @test all(tracer ∈ keys(model.tracers) for tracer in required_tracers)
    @test :PAR ∈ keys(model.auxiliary_fields) && isa(model.auxiliary_fields.PAR, TwoBandPhotosyntheticallyActiveRatiationField)

    # checks model works with zero values
    time_step!(model, 1.0)

    # and that they all return zero
    @test all([all(values .== 0) for values in values(model.tracers)]) 
    
    return nothing
end

for arch in (CPU(), )
    grid = RectilinearGrid(arch; size=(3, 3, 6), extent=(1, 1, 2))
    for carbonates = (false, true), oxygen = (false, true), sinking = (false, true), open_bottom = (false, true)
        if !(sinking && open_bottom) # no sinking is the same with and without open bottom
            @info "Testing on $arch with carbonates $(carbonates ? :✅ : :❌), oxygen $(oxygen ? :✅ : :❌), sinking $(sinking ? :✅ : :❌), open bottom $(open_bottom ? :✅ : :❌))"
            @testset "$arch, $carbonates, $oxygen, $sinking, $open_bottom" begin
                test_LOBSTER(grid, carbonates, oxygen, sinking, open_bottom)
            end
        end
    end
end
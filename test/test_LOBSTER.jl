using Test
using OceanBioME: LOBSTER
using Oceananigans

function test_LOBSTER(grid, carbonates, oxygen, variable_redfield, sinking, open_bottom)
    PAR = CenterField(grid)

    if sinking
        model = NonhydrostaticModel(;grid,
                                     biogeochemistry = LOBSTER(;grid, carbonates, oxygen, variable_redfield, open_bottom),
                                     auxiliary_fields = (; PAR))
    else
        model = NonhydrostaticModel(;grid,
                                     biogeochemistry = LOBSTER(;grid, carbonates, oxygen, variable_redfield, sinking_velocities = NamedTuple()),
                                     auxiliary_fields = (; PAR))
    end

    # correct tracers and auxiliary fields have been setup, and order has not changed
    required_tracers = (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM)
    if carbonates
        required_tracers = (required_tracers..., :DIC, :ALK)
    end
    if oxygen
        required_tracers = (required_tracers..., :OXY)
    end
    if variable_redfield
        required_tracers = (required_tracers..., :Dᶜ, :DDᶜ, :DOMᶜ)
    end

    @test Oceananigans.Biogeochemistry.required_biogeochemical_tracers(model.biogeochemistry) == required_tracers
    @test all(tracer ∈ keys(model.tracers) for tracer in required_tracers)
    @test :PAR ∈ keys(model.auxiliary_fields)

    # checks model works with zero values
    time_step!(model, 1.0)

    # and that they all return zero
    @test all([all(values .== 0) for values in values(model.tracers)]) 

    # mass conservation
    model.tracers.NO₃ .= rand()
    model.tracers.NH₄ .= rand()
    model.tracers.P .= rand()
    model.tracers.Z .= rand()
    model.tracers.D .= rand()
    model.tracers.DD .= rand()
    model.tracers.DOM .= rand()

    ΣN₀ = sum(model.tracers.NO₃) + sum(model.tracers.NH₄) + sum(model.tracers.P) + sum(model.tracers.Z) + sum(model.tracers.D) + sum(model.tracers.DD) + sum(model.tracers.DOM)

    time_step!(model, 1.0)

    ΣN₁ = sum(model.tracers.NO₃) + sum(model.tracers.NH₄) + sum(model.tracers.P) + sum(model.tracers.Z) + sum(model.tracers.D) + sum(model.tracers.DD) + sum(model.tracers.DOM)
    
    @test ΣN₀ ≈ ΣN₁ # guess this should actually fail with a high enough accuracy when sinking is on with an open bottom

    return nothing
end

for arch in (CPU(), )
    grid = RectilinearGrid(arch; size=(3, 3, 6), extent=(1, 1, 2))
    for carbonates = (false, true), oxygen = (false, true), variable_redfield = (true, false), sinking = (false, true), open_bottom = (false, true)
        if !(sinking && open_bottom) # no sinking is the same with and without open bottom
            @info "Testing on $arch with carbonates $(carbonates ? :✅ : :❌), oxygen $(oxygen ? :✅ : :❌), variable redfield $(variable_redfield ? :✅ : :❌), sinking $(sinking ? :✅ : :❌), open bottom $(open_bottom ? :✅ : :❌))"
            @testset "$arch, $carbonates, $oxygen, $variable_redfield, $sinking, $open_bottom" begin
                test_LOBSTER(grid, carbonates, oxygen, variable_redfield, sinking, open_bottom)
            end
        end
    end
end
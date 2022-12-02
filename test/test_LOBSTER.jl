using Test
using OceanBioME: LOBSTER
using Oceananigans

ΣN(model, variable_redfield) = variable_redfield ? (
    sum(model.tracers.NO₃) + sum(model.tracers.NH₄) + sum(model.tracers.P) + sum(model.tracers.Z) + sum(model.tracers.sPON) + sum(model.tracers.bPON) + sum(model.tracers.DON)) : (
    sum(model.tracers.NO₃) + sum(model.tracers.NH₄) + sum(model.tracers.P) + sum(model.tracers.Z) + sum(model.tracers.sPOM) + sum(model.tracers.bPOM) + sum(model.tracers.DOM))

function ΣC(model, carbonates, variable_redfield)
    # *will only be conserved if carbonates on*
    if variable_redfield
        OC = sum(model.tracers.sPOC .+ model.tracers.bPOC .+ model.tracers.DOC)
    else
        OC = sum(model.tracers.sPOM .+ model.tracers.bPOM .+ model.tracers.DOM) * model.biogeochemistry.disolved_organic_redfield
    end
    
    if carbonates
        IC = sum(model.tracers.DIC)
    else
        IC = 0.0
    end

    LC = (model.tracers.P .+ model.tracers.Z) * model.biogeochemistry.phytoplankton_redfield

    return OC + IC + LC
end

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
    required_tracers = variable_redfield ? (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON) : (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)
    if carbonates
        required_tracers = (required_tracers..., :DIC, :Alk)
    end
    if oxygen
        required_tracers = (required_tracers..., :O₂)
    end
    if variable_redfield
        required_tracers = (required_tracers..., :sPOC, :bPOC, :DOC)
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
    if variable_redfield
        model.tracers.sPON .= rand()
        model.tracers.bPON .= rand()
        model.tracers.DON .= rand()
    else
        model.tracers.sPOM .= rand()
        model.tracers.bPOM .= rand()
        model.tracers.DOM .= rand()
    end

    ΣN₀ = ΣN(model, variable_redfield)

    ΣC₀ = ΣC(model, carbonates, variable_redfield)

    time_step!(model, 1.0)

    ΣN₁ = ΣN(model, variable_redfield)
    
    @test ΣN₀ ≈ ΣN₁ # guess this should actually fail with a high enough accuracy when sinking is on with an open bottom

    if carbonates
        ΣC₁ = ΣC(model, carbonates, variable_redfield)

        @test ΣC₀ ≈ ΣC₁
    end

    return nothing
end

for arch in (CPU(), )
    grid = RectilinearGrid(arch; size=(3, 3, 6), extent=(1, 1, 2))
    for carbonates = (false, true), oxygen = (false, true), variable_redfield = (false, true), sinking = (false, true), open_bottom = (false, true)
        if !(sinking && open_bottom) # no sinking is the same with and without open bottom
            @info "Testing on $arch with carbonates $(carbonates ? :✅ : :❌), oxygen $(oxygen ? :✅ : :❌), variable redfield $(variable_redfield ? :✅ : :❌), sinking $(sinking ? :✅ : :❌), open bottom $(open_bottom ? :✅ : :❌))"
            @testset "$arch, $carbonates, $oxygen, $variable_redfield, $sinking, $open_bottom" begin
                test_LOBSTER(grid, carbonates, oxygen, variable_redfield, sinking, open_bottom)
            end
        end
    end
end
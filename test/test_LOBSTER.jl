using Test
using OceanBioME: LOBSTER, conserved_tracers, redfield
using Oceananigans, CUDA

ΣN(model, biogeochemistry) = sum([sum(model.tracers[tracer_name]) for tracer_name in conserved_tracers(biogeochemistry)])

function ΣC(model, carbonates, biogeochemistry)
    OC = sum([sum(model.tracers[tracer_name] * redfield(Val(tracer_name), biogeochemistry, model.tracers)) for tracer_name in conserved_tracers(biogeochemistry)])

    if carbonates
        IC = sum(model.tracers.DIC)
    else
        IC = 0.0
    end

    return OC + IC
end

ΣGⁿ(model, biogeochemistry) = sum([sum(model.timestepper.Gⁿ[tracer_name]) for tracer_name in conserved_tracers(biogeochemistry)])


function ΣGᶜ(model, biogeochemistry)
    # *will only be conserved if carbonates on*
    # Leaving this as is for now since as the `redfield` is based on value not tendency then it returns the wrong results if
    # the redfield is rapidly changing
    OC = sum(model.timestepper.Gⁿ.sPOC .+ model.timestepper.Gⁿ.bPOC .+ model.timestepper.Gⁿ.DOC)

    IC = sum(model.timestepper.Gⁿ.DIC)

    LC = sum(model.timestepper.Gⁿ.P * (1 + biogeochemistry.organic_carbon_calcate_ratio) .+ model.timestepper.Gⁿ.Z) * biogeochemistry.phytoplankton_redfield 

    return OC + IC + LC
end

function test_LOBSTER(grid, carbonates, oxygen, variable_redfield, sinking, open_bottom, n_timesteps)
    if sinking
        model = NonhydrostaticModel(;grid,
                                     biogeochemistry = LOBSTER(;grid, carbonates, oxygen, variable_redfield, open_bottom))
    else
        model = NonhydrostaticModel(;grid,
                                     biogeochemistry = LOBSTER(;grid, carbonates, oxygen, variable_redfield, sinking_speeds = NamedTuple()))
    end

    # correct tracers and auxiliary fields have been setup, and order has not changed
    required_tracers = conserved_tracers(model.biogeochemistry)
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

    # checks model works with zero values
    time_step!(model, 1.0)

    # mass conservation
    CUDA.@allowscalar begin 
        # and that they all return zero
        @test all([all(values .== 0) for values in values(model.tracers)]) 

        model.tracers.NO₃ .= 10 * rand()
        model.tracers.NH₄ .= rand()
        model.tracers.P .= rand()
        model.tracers.Z .= rand()
        if variable_redfield
            model.tracers.sPON .= rand()
            model.tracers.bPON .= rand()
            model.tracers.DON .= rand()

            model.tracers.sPOC .= model.tracers.sPON * redfield(Val(:sPOM), model.biogeochemistry)
            model.tracers.bPOC .= model.tracers.bPON * redfield(Val(:bPOM), model.biogeochemistry)
            model.tracers.DOC .= model.tracers.DON * redfield(Val(:DOM), model.biogeochemistry)
        else
            model.tracers.sPOM .= rand()
            model.tracers.bPOM .= rand()
            model.tracers.DOM .= rand()
        end

        if carbonates
            model.tracers.DIC .= 2000 * rand()
            model.tracers.Alk .= 2000 * rand()
        end
    end
    
    ΣN₀ = CUDA.@allowscalar ΣN(model, model.biogeochemistry)

    ΣC₀ = CUDA.@allowscalar ΣC(model, carbonates, model.biogeochemistry)
    
    for _ in 1:n_timesteps
        time_step!(model, 1.0)
    end

    ΣN₁ = CUDA.@allowscalar ΣN(model, model.biogeochemistry)
    
    CUDA.@allowscalar if !(sinking && open_bottom) #when we have open bottom sinking we won't conserve anything
        @test ΣN₀ ≈ ΣN₁
        @test ΣGⁿ(model, model.biogeochemistry) ≈ 0.0 atol = 1e-15 # rtol=sqrt(eps) so is usually much larger than even this

        if (carbonates && variable_redfield)
            ΣC₁ = ΣC(model, carbonates, model.biogeochemistry)
            @test ΣC₀ ≈ ΣC₁

            @test ΣGᶜ(model, model.biogeochemistry.underlying_biogeochemistry) ≈ 0.0 atol = 1e-15
        end
    end
    return model
end

n_timesteps = 100

grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

for open_bottom = (false, true), sinking = (false, true), variable_redfield = (false, true), oxygen = (false, true), carbonates = (false, true)
    if !(sinking && open_bottom) # no sinking is the same with and without open bottom
        @info "Testing on $(typeof(architecture)) with carbonates $(carbonates ? :✅ : :❌), oxygen $(oxygen ? :✅ : :❌), variable redfield $(variable_redfield ? :✅ : :❌), sinking $(sinking ? :✅ : :❌), open bottom $(open_bottom ? :✅ : :❌))"
        @testset "$architecture, $carbonates, $oxygen, $variable_redfield, $sinking, $open_bottom" begin
            test_LOBSTER(grid, carbonates, oxygen, variable_redfield, sinking, open_bottom, n_timesteps)
        end
    end
end
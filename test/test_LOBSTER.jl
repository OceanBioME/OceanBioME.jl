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
        OC = sum(model.tracers.sPOM .+ model.tracers.bPOM .+ model.tracers.DOM) * model.biogeochemistry. organic_redfield
    end
    
    if carbonates
        IC = sum(model.tracers.DIC)
    else
        IC = 0.0
    end

    LC = sum(model.tracers.P * (1 + model.biogeochemistry.organic_carbon_calcate_ratio) .+ model.tracers.Z) * model.biogeochemistry.phytoplankton_redfield 

    return OC + IC + LC
end

ΣGⁿ(model, variable_redfield) = variable_redfield ? (
    sum(model.timestepper.Gⁿ.NO₃) + sum(model.timestepper.Gⁿ.NH₄) + sum(model.timestepper.Gⁿ.P) + sum(model.timestepper.Gⁿ.Z) + sum(model.timestepper.Gⁿ.sPON) + sum(model.timestepper.Gⁿ.bPON) + sum(model.timestepper.Gⁿ.DON)) : (
    sum(model.timestepper.Gⁿ.NO₃) + sum(model.timestepper.Gⁿ.NH₄) + sum(model.timestepper.Gⁿ.P) + sum(model.timestepper.Gⁿ.Z) + sum(model.timestepper.Gⁿ.sPOM) + sum(model.timestepper.Gⁿ.bPOM) + sum(model.timestepper.Gⁿ.DOM))


function ΣGᶜ(model, carbonates, variable_redfield)
    # *will only be conserved if carbonates on*
    if variable_redfield
        OC = sum(model.timestepper.Gⁿ.sPOC .+ model.timestepper.Gⁿ.bPOC .+ model.timestepper.Gⁿ.DOC)
    else
        OC = sum(model.timestepper.Gⁿ.sPOM .+ model.timestepper.Gⁿ.bPOM .+ model.timestepper.Gⁿ.DOM) * model.biogeochemistry. organic_redfield
    end
    
    if carbonates
        IC = sum(model.timestepper.Gⁿ.DIC)
    else
        IC = 0.0
    end

    LC = sum(model.timestepper.Gⁿ.P .+ model.timestepper.Gⁿ.Z) * model.biogeochemistry.phytoplankton_redfield 

    return OC + IC + LC
end
function test_LOBSTER(grid, carbonates, oxygen, variable_redfield, sinking, open_bottom, n_timesteps)
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
    model.tracers.NO₃ .= 10 * rand()
    model.tracers.NH₄ .= rand()
    model.tracers.P .= rand()
    model.tracers.Z .= rand()
    if variable_redfield
        model.tracers.sPON .= rand()
        model.tracers.bPON .= rand()
        model.tracers.DON .= rand()

        model.tracers.sPOC .= model.tracers.sPON * model.biogeochemistry.phytoplankton_redfield
        model.tracers.bPOC .= model.tracers.bPON * model.biogeochemistry.phytoplankton_redfield
        model.tracers.DOC .= model.tracers.DON * model.biogeochemistry.phytoplankton_redfield
    else
        model.tracers.sPOM .= rand()
        model.tracers.bPOM .= rand()
        model.tracers.DOM .= rand()
    end

    if carbonates
        model.tracers.DIC .= 2000 * rand()
        model.tracers.Alk .= 200 * rand()
    end

    ΣN₀ = ΣN(model, variable_redfield)

    ΣC₀ = ΣC(model, carbonates, variable_redfield)
    
    for i in 1:n_timesteps
        time_step!(model, 1.0)
    end

    ΣN₁ = ΣN(model, variable_redfield)
    
    if !(sinking && open_bottom) #when we have open bottom sinking we won't conserve anything
        @test ΣN₀ ≈ ΣN₁
        @test ΣGⁿ(model, variable_redfield) ≈ 0.0 atol = 1e-15 # rtol=sqrt(eps) so is usually much larger than even this

        if carbonates
            ΣC₁ = ΣC(model, carbonates, variable_redfield)
            @test ΣC₀ ≈ ΣC₁# atol = 0.0001 # when we convert to and from 

            @test ΣGᶜ(model, carbonates, variable_redfield) ≈ 0.0 atol = 1e-15
        end
    end
    return model
end

n_timesteps = 1

for arch in (CPU(), )
    grid = RectilinearGrid(arch; size=(1, 1, 1), extent=(1, 1, 2))
    for open_bottom = (false, true), sinking = (false, true), variable_redfield = (false, true), oxygen = (false, true), carbonates = (false, true)
        if !(sinking && open_bottom) # no sinking is the same with and without open bottom
            @info "Testing on $arch with carbonates $(carbonates ? :✅ : :❌), oxygen $(oxygen ? :✅ : :❌), variable redfield $(variable_redfield ? :✅ : :❌), sinking $(sinking ? :✅ : :❌), open bottom $(open_bottom ? :✅ : :❌))"
            @testset "$arch, $carbonates, $oxygen, $variable_redfield, $sinking, $open_bottom" begin
                test_LOBSTER(grid, carbonates, oxygen, variable_redfield, sinking, open_bottom, n_timesteps)
            end
        end
    end
end
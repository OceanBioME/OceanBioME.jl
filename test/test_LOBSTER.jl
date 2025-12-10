#include("dependencies_for_runtests.jl")

using OceanBioME: conserved_tracers
using Oceananigans, CUDA, Random

Random.seed!(42)

ΣN(model, biogeochemistry) = sum([sum(model.tracers[tracer_name]) for tracer_name in conserved_tracers(biogeochemistry, true).nitrogen])

function ΣC(model, carbonates, biogeochemistry)
    conserved = conserved_tracers(biogeochemistry, true)
    
    if !(:carbon in keys(conserved))
        return NaN
    end

    conserved_tracer_names = conserved.carbon.tracers
    conserved_tracer_scalefactors = conserved.carbon.scalefactors

    return sum([sum(model.tracers[tracer_name]) * conserved_tracer_scalefactors[n] for (n, tracer_name) in enumerate(conserved_tracer_names)])
end

ΣGⁿ(model, biogeochemistry) = sum([sum(model.timestepper.Gⁿ[tracer_name]) for tracer_name in conserved_tracers(biogeochemistry, true).nitrogen])

function ΣGᶜ(model, biogeochemistry)
    conserved = conserved_tracers(biogeochemistry, true)
    
    if !(:carbon in keys(conserved))
        return NaN
    end

    conserved_tracer_names = conserved.carbon.tracers
    conserved_tracer_scalefactors = conserved.carbon.scalefactors

    return sum([sum(model.timestepper.Gⁿ[tracer_name]) * conserved_tracer_scalefactors[n] for (n, tracer_name) in enumerate(conserved_tracer_names)])
end

function test_LOBSTER(grid, nutrients, carbonate_system, oxygen, detritus, sinking, open_bottom, n_timesteps)
    model = NonhydrostaticModel(; grid,
                                  biogeochemistry = LOBSTER(; grid, nutrients, carbonate_system, oxygen, detritus))

    # correct tracers and auxiliary fields have been setup, and order has not changed
    required_tracers = (:NO₃, :NH₄, :P, :Z)
    if nutrients isa NitrateAmmoniaIron
        required_tracers = (required_tracers..., :Fe)
    end
    if !isnothing(carbonate_system)
        required_tracers = (required_tracers..., :DIC, :Alk)
    end
    if !isnothing(oxygen)
        required_tracers = (required_tracers..., :O₂)
    end
    if detritus isa VariableRedfieldDetritus
        required_tracers = (required_tracers..., :sPON, :bPON, :DON, :sPOC, :bPOC, :DOC)
    else
        required_tracers = (required_tracers..., :sPOM, :bPOM, :DOM)
    end

    @test all(tracer ∈ Oceananigans.Biogeochemistry.required_biogeochemical_tracers(model.biogeochemistry) for tracer in required_tracers)
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
        if detritus isa VariableRedfieldDetritus
            model.tracers.sPON .= rand()
            model.tracers.bPON .= rand()
            model.tracers.DON .= rand()

            model.tracers.sPOC .= model.tracers.sPON * 6.56
            model.tracers.bPOC .= model.tracers.bPON * 6.56
            model.tracers.DOC .= model.tracers.DON * 6.56
        else
            model.tracers.sPOM .= rand()
            model.tracers.bPOM .= rand()
            model.tracers.DOM .= rand()
        end

        if !isnothing(carbonate_system)
            model.tracers.DIC .= 2000 * rand()
            model.tracers.Alk .= 2000 * rand()
        end
    end
    
    ΣN₀ = CUDA.@allowscalar ΣN(model, model.biogeochemistry)

    ΣC₀ = CUDA.@allowscalar ΣC(model, carbonate_system, model.biogeochemistry)

    for _ in 1:n_timesteps
        time_step!(model, 1.0)
    end

    ΣN₁ = CUDA.@allowscalar ΣN(model, model.biogeochemistry)
    
    CUDA.@allowscalar if !(sinking && open_bottom) #when we have open bottom sinking we won't conserve anything
        @test ΣN₀ ≈ ΣN₁
        @test ΣGⁿ(model, model.biogeochemistry) ≈ 0.0 atol = 1e-15 # rtol=sqrt(eps) so is usually much larger than even this

        if !isnothing(carbonate_system)
            ΣC₁ = ΣC(model, carbonate_system, model.biogeochemistry)
            @test ΣC₀ ≈ ΣC₁
            @test ΣGᶜ(model, model.biogeochemistry.underlying_biogeochemistry) ≈ 0.0 atol = 1e-15
        end
    end
    return model
end

n_timesteps = 100

grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

for open_bottom = (false, true), 
    sinking = (false, true), 
    detritus = (TwoParticleAndDissolved, VariableRedfieldDetritus),
    oxygen = (nothing, Oxygen()), 
    carbonate_system = (nothing, CarbonateSystem()),
    nutrients = (NitrateAmmonia(), NitrateAmmoniaIron())

    detritus = detritus(grid; 
                        small_particle_sinking_speed = sinking ? 3.47e-5 : 0.0, 
                        large_particle_sinking_speed = sinking ? 200/day : 0.0,
                        open_bottom)

    if !(sinking && open_bottom) # no sinking is the same with and without open bottom
        @testset "LOBSTER ($(nameof(typeof(nutrients)))/$(nameof(typeof(carbonate_system)))/$(nameof(typeof(oxygen)))/$(nameof(typeof(detritus)))/$sinking/$open_bottom" begin
            test_LOBSTER(grid, nutrients, carbonate_system, oxygen, detritus, sinking, open_bottom, n_timesteps)
        end
    end
end

@testset "Float32 LOBSTER" begin
    grid = RectilinearGrid(architecture, Float32; size=(3, 3, 10), extent=(10, 10, 200))
    bgc = LOBSTER(; grid)

    par = bgc.underlying_biogeochemistry.light_attenuation
    @test par.water_red_attenuation isa Float32
    @test par.water_blue_attenuation isa Float32
    @test par.chlorophyll_red_attenuation isa Float32
    @test par.chlorophyll_blue_attenuation isa Float32
    @test par.chlorophyll_red_exponent isa Float32
    @test par.chlorophyll_blue_exponent isa Float32
    @test par.pigment_ratio isa Float32
    @test par.phytoplankton_chlorophyll_ratio isa Float32
end
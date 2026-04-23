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

n_timesteps = 100

grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

instantiate_detritus(::Val{detritus}, grid, sinking) where detritus = 
    detritus(grid;
             small_particle_sinking_speed = sinking ? 3.47e-5 : 0.0,
             large_particle_sinking_speed = sinking ? 200/day : 0.0)

instantiate_detritus(::Val{Detritus}, grid, sinking) = 
    Detritus(grid;
             sinking_speed = sinking ? 0.1213 / day : 0.0)

test_models = []

light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))

# construct the possible configurations
for sinking = (false, true),
    detritus = (TwoParticleAndDissolved, VariableRedfieldDetritus, Detritus),
    oxygen = (nothing, Oxygen()),
    carbonate_system = (nothing, CarbonateSystem()),
    nutrients = (NitrateAmmonia(), NitrateAmmoniaIron(), Nutrient())

    detritus = instantiate_detritus(Val(detritus), grid, sinking)

    model = NonhydrostaticModel(grid;
                                biogeochemistry = NutrientsBiologyDetritus(grid; nutrients, 
                                                                                biology = PhytoZoo(),
                                                                                carbonate_system, 
                                                                                oxygen, 
                                                                                detritus,
                                                                                light_attenuation))

    required_tracers = (:P, :Z)
    initial_values = (rand(), rand())

    if nutrients isa NitrateAmmonia
        required_tracers = (required_tracers..., :NO₃, :NH₄)
        initial_values = (initial_values..., 10*rand(), rand())
    elseif nutrients isa NitrateAmmoniaIron
        required_tracers = (required_tracers..., :NO₃, :NH₄, :Fe)
        initial_values = (initial_values..., 10*rand(), rand(), 1e-3*rand())
    elseif nutrients isa Nutrient
        required_tracers = (required_tracers..., :N)
        initial_values = (initial_values..., 10*rand())
    end

    if !isnothing(carbonate_system)
        required_tracers = (required_tracers..., :DIC, :Alk)
        initial_values = (initial_values..., 2000*rand(), 2000*rand())
    end

    if !isnothing(oxygen)
        required_tracers = (required_tracers..., :O₂)
        initial_values = (initial_values..., 300*rand())
    end

    if detritus isa VariableRedfieldDetritus
        required_tracers = (required_tracers..., :sPON, :bPON, :DON, :sPOC, :bPOC, :DOC)
        sPON, bPON, DON = rand(), rand(), rand()
        initial_values = (initial_values..., sPON, bPON, DON, 6.56*sPON, 6.56*bPON, 6.56*DON)
    elseif detritus isa TwoParticleAndDissolved
        required_tracers = (required_tracers..., :sPOM, :bPOM, :DOM)
        initial_values = (initial_values..., rand(), rand(), rand())
    else
        required_tracers = (required_tracers..., :D)
        initial_values = (initial_values..., rand())
    end

    @info "Constructed $(keys(model.tracers))"

    push!(test_models, (; required_tracers, initial_values, model, sinking))
end

# since compile time dominates tests we can test different conservations separatly
# at little extra cost because once its been stepped its ~instantanous to step again

@testset "LOBSTER setup and dry step" (for (required_tracers, initial_values, model, sinking) in test_models
    tracers_are_correct = all(tracer ∈ Oceananigans.Biogeochemistry.required_biogeochemical_tracers(model.biogeochemistry) for tracer in required_tracers)
    tracers_are_materialised =  all(tracer ∈ keys(model.tracers) for tracer in required_tracers)

    time_step!(model, 1.0) # long wait for compile here...

    # CUDA.@allowscalar maybe
    everything_stayed_zero = CUDA.@allowscalar all([all(values .== 0) for values in values(model.tracers)])

    @test tracers_are_correct & tracers_are_materialised & everything_stayed_zero
    @info "Stepped $(keys(model.tracers))"
end; )

@testset "LOBSTER nitrogen conservation" (for (required_tracers, initial_values, model, sinking) in test_models
    if !sinking
        @info "Testing $(keys(model.tracers))"
        set!(model; NamedTuple{required_tracers}(initial_values)...)
        ΣN₀ = CUDA.@allowscalar ΣN(model, model.biogeochemistry)
        for _ in 1:n_timesteps
            time_step!(model, 1.0)
        end
        ΣN₁ = CUDA.@allowscalar ΣN(model, model.biogeochemistry)
        @test ΣN₀ ≈ ΣN₁
    end
end; )

@testset "LOBSTER carbon conservation" (for (required_tracers, initial_values, model, sinking) in test_models
    carbonate_system = model.biogeochemistry.underlying_biogeochemistry.carbonate_system
    if !(isnothing(carbonate_system)|sinking)
        @info "Testing $(keys(model.tracers))"
        set!(model; NamedTuple{required_tracers}(initial_values)...)
        ΣC₀ = CUDA.@allowscalar ΣC(model, carbonate_system, model.biogeochemistry)
        for _ in 1:n_timesteps
            time_step!(model, 1.0)
        end
        ΣC₁ = ΣC(model, carbonate_system, model.biogeochemistry)
        @test ΣC₀ ≈ ΣC₁
    end
end; )

@testset "BND constructors" begin
    grid = RectilinearGrid(architecture; size=(1, 1, 1), extent=(1, 1, 2))

    lobster = LOBSTER(grid)
    npzd = NPZD(grid)

    @test lobster isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsBiologyDetritus}
    @test npzd isa OceanBioME.DiscreteBiogeochemistry{<:NutrientsBiologyDetritus}

    @test lobster.underlying_biogeochemistry isa NutrientsBiologyDetritus{<:NitrateAmmonia, <:PhytoZoo, <:TwoParticleAndDissolved}
    @test npzd.underlying_biogeochemistry isa NutrientsBiologyDetritus{<:Nutrient, <:PhytoZoo, <:Detritus}
end


@testset "Float32 LOBSTER" begin
    grid = RectilinearGrid(architecture, Float32; size=(3, 3, 10), extent=(10, 10, 200))
    bgc = LOBSTER(grid)

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

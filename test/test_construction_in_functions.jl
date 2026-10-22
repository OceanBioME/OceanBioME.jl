include("dependencies_for_runtests.jl")

using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity

# These components used to define their tendency methods with `eval` when constructed, which the
# function that called the constructor cannot see (Julia's world-age rule): a model built and run
# inside a function silently lost tendencies and sinking, or errored (#425). So every model here is
# built, stepped and inspected inside a function. The tests only catch that regression if they run
# before any other test constructs the same components at top level (hence their place in runtests.jl).

grid = RectilinearGrid(architecture; size = (1, 1, 8), extent = (1, 1, 80))

light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(ConstantField(100))

total(model, names) = sum(sum(interior(model.tracers[name])) for name in names)

top_cell(model, name) = CUDA.@allowscalar model.tracers[name][1, 1, model.grid.Nz]

# builds a LOBSTER column with the given detritus class names, steps it, and reports what the model
# did, all inside this function
function run_lobster_inside_function(grid, dissolved, particulate)
    detritus = DissolvedParticulate(grid, dissolved, particulate; open_bottom = false)
    plankton = PhytoZoo(; edible_detritus_name = first(particulate)) # `PhytoZoo` grazes one class by name
    biogeochemistry = LOBSTER(grid; detritus, plankton, light_attenuation)
    model = NonhydrostaticModel(grid; biogeochemistry)

    # with no plankton the detritus only remineralises (into NH₄, via the dissolved class) and sinks,
    # and with a closed bottom that conserves nitrogen
    detritus_names = (dissolved, particulate...)
    set!(model; NO₃ = 5, NH₄ = 0.1, NamedTuple{detritus_names}(ntuple(_ -> 1, length(detritus_names)))...)

    nitrogen_names = (:NO₃, :NH₄, :P, :Z, detritus_names...)
    initial_nitrogen = total(model, nitrogen_names)

    for _ in 1:10
        time_step!(model, 60)
    end

    return (sinks = !isnothing(biogeochemical_drift_velocity(model.biogeochemistry, Val(first(particulate)))),
            particulate_decayed = top_cell(model, first(particulate)) < 1,
            nitrogen_conserved = isapprox(total(model, nitrogen_names), initial_nitrogen; rtol = 1e-12))
end

@testset "DissolvedParticulate constructed and used inside a function" begin
    # the default LOBSTER class names, and names no other test uses
    for (dissolved, particulate) in ((:DOM, (:sPOM, :bPOM)), (:DOMᵢ, (:sPOMᵢ, :bPOMᵢ)))
        results = run_lobster_inside_function(grid, dissolved, particulate)

        @test results.sinks
        @test results.particulate_decayed
        @test results.nitrogen_conserved
    end
end

implicit_replicates(grid) = CarbonateSystem(2)
explicit_replicates(grid) = ExplicitCalciumCarbonate(grid; replicates = 2)
labelled_implicit_replicates(grid) = CarbonateSystem((:a, :b))
labelled_explicit_replicates(grid) = ExplicitCalciumCarbonate(grid; replicates = (:a, :b))

# builds an NPZD column with two identically initialised replicates of the carbonate tracers, each with
# its own air-sea CO₂ exchange, steps it, and reports whether the first replicate responded to the
# primary production and whether the replicates are still identical, all inside this function
function run_replicated_carbonate_inside_function(grid, build_inorganic_carbon)
    inorganic_carbon = build_inorganic_carbon(grid)
    biogeochemistry = NPZD(grid; inorganic_carbon, light_attenuation)
    boundary_conditions = carbon_dioxide_gas_exchange_boundary_conditions(biogeochemistry)

    # the carbon chemistry behind `ExplicitCalciumCarbonate` and the gas exchange needs salinity as well as temperature
    model = NonhydrostaticModel(grid; biogeochemistry, boundary_conditions, auxiliary_fields = (; S = ConstantField(35.0)))

    a, b = carbonate_replicates(biogeochemistry)

    set!(model; N = 5, P = 1, Z = 0.5, D = 1, T = 15, (a.DIC => 2100, a.Alk => 2300, b.DIC => 2100, b.Alk => 2300)...)

    for _ in 1:10
        time_step!(model, 60)
    end

    identical(name) = Array(interior(fields(model)[a[name]])) == Array(interior(fields(model)[b[name]]))

    return (replicate_changed = top_cell(model, a.DIC) != 2100,
            replicates_identical = all(identical, filter(!=(:label), keys(a))),
            calcium_carbonate_sinks = haskey(a, :CaCO₃) && !isnothing(biogeochemical_drift_velocity(model.biogeochemistry, Val(a.CaCO₃))))
end

@testset "Replicated carbonate tracers constructed and used inside a function" begin
    for build_inorganic_carbon in (implicit_replicates, labelled_implicit_replicates)
        implicit = run_replicated_carbonate_inside_function(grid, build_inorganic_carbon)

        @test implicit.replicate_changed
        @test implicit.replicates_identical
    end

    for build_inorganic_carbon in (explicit_replicates, labelled_explicit_replicates)
        explicit = run_replicated_carbonate_inside_function(grid, build_inorganic_carbon)

        @test explicit.replicate_changed
        @test explicit.replicates_identical
        @test explicit.calcium_carbonate_sinks
    end
end

sediment_grid = RectilinearGrid(architecture; size = (3, 3, 10), extent = (1, 1, 100), halo = (3, 3, 3))

# builds an NPZD model with an instant-remineralisation sediment, steps it, and reports whether anything
# reached the sediment, all inside this function
function run_sediment_inside_function(grid)
    biogeochemistry = NPZD(grid; sediment = InstantRemineralisationSediment(grid), light_attenuation)
    model = NonhydrostaticModel(grid; biogeochemistry)

    set!(model; N = 1, P = 1, Z = 1, D = 1)

    for _ in 1:10
        time_step!(model, 60)
    end

    return (; buried = sum(interior(biogeochemistry.sediment.fields.storage)) > 0)
end

@testset "InstantRemineralisation sediment constructed and used inside a function" begin
    results = run_sediment_inside_function(sediment_grid)

    @test results.buried
end

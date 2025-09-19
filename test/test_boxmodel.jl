#include("dependencies_for_runtests.jl")

using OceanBioME.BoxModels: boxmodel_xyz

using Oceananigans: AbstractGrid, AbstractModel, AbstractOutputWriter
using Oceananigans.Fields: FunctionField
using Oceananigans.Grids: nodes

PAR(t) = (60 * (1 - cos((t + 15days) * 2π / 365days)) * (1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days)^2))) + 2) * exp(-2)

defaults = (NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01)

function simple_box_model(; x = nothing, y = nothing, z = nothing)
    grid = BoxModelGrid(; x, y, z)

    clock = Clock(time = zero(grid))

    light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(FunctionField{Center, Center, Center}(PAR, grid; clock))

    biogeochemistry = LOBSTER(; grid, light_attenuation)

    model = BoxModel(; grid, biogeochemistry, clock)

    return model
end

@testset "Construct `BoxModel`" begin
    for z in (nothing, -100)
        model = simple_box_model(; z)

        @test model.grid isa AbstractGrid
        @test model isa AbstractModel
        @test set!(model; defaults...) isa Nothing
        @test boxmodel_xyz(nodes(model.grid, Center(), Center(), Center()), model.grid) == (0, 0, ifelse(isnothing(z), 0, z))
    end
end

@testset "Timestep `BoxModel`" begin
    model = simple_box_model()

    set!(model; defaults...) 

    for n in 1:10
        time_step!(model, 20minutes)
    end

    # values are being updated
    for (name, field) in pairs(model.fields)
        default = name in keys(defaults) ? defaults[name] : 0

        @test field[1, 1, 1] != default
    end
end

@testset "`BoxModel` simulation" begin
    model = simple_box_model()

    set!(model; defaults...) 

    simulation = Simulation(model; Δt = 20minutes, stop_iteration = 1000, verbose = false)

    fname = "box_model_test.jld2"

    fast_output = SpeedyOutput(fname)

    @test_broken fast_output isa AbstractOutputWriter

    simulation.callbacks[:output] = Callback(fast_output, IterationInterval(20))

    @test run!(simulation) isa Nothing

    results = load_output(fast_output)

    # outputting as expected
    @test length(keys(results)) == length(model.fields) + 1
    @test length(results.t) == 1000 / 20 + 1

    # something is happening and we're not just writing the same value 
    @test all([result[1] != result[end] for result in values(results)])
end
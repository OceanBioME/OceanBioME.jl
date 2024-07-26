#include("dependencies_for_runtests.jl")

using Oceananigans: AbstractGrid, AbstractModel, AbstractOutputWriter
using Oceananigans.Fields: FunctionField

PAR(t) = (60 * (1 - cos((t + 15days) * 2π / 365days)) * (1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days)^2))) + 2) * exp(-2)

defaults = (NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01)

function simple_box_model()
    grid = BoxModelGrid()
    clock = Clock(time = zero(grid))

    light_attenuation_model = PrescribedPhotosyntheticallyActiveRadiation(FunctionField{Center, Center, Center}(PAR, grid; clock))

    biogeochemistry = LOBSTER(; grid, light_attenuation_model)

    model = BoxModel(; grid, biogeochemistry, clock)

    return model
end

@testset "Construct `BoxModel`" begin
    model = simple_box_model()

    @test grid isa AbstractGrid
    @test model isa AbstractModel
    @test set!(model; defaults...) isa Nothing
end

@testset "Timestep `BoxModel`" begin
    model = simple_box_model()

    set!(model; defaults...) 

    for n in 1:10
        time_step!(model, 20minutes)
    end

    # values are being updated
    for (name, field) in enumerate(model.fields)
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
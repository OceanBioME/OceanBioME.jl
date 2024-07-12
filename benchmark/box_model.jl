# # Calibrating a biogeochemical model with `EnsembleKalmanProcesses`
#
# In this example we calibrate some of the parameters for the [NPZD](@ref NPZD) model
# in a simple box model setup using a data assimilation package [EnsembleKalmanProcesses](https://github.com/CliMA/EnsembleKalmanProcesses.jl).
# First we setup the model and generate synthetic data with "true" parameters. We then
# define priors and setup an EKP to solve.
#
# While this is a very simple situation it illustrates the ease of integration with 
# data assimilation tools. Examples given in the EnsembleKalmanProcesses docs illustrate
# how the package can be used to solve more complex forward models.

# ## Install dependencies
# First we ensure we have the required dependencies installed
# ```julia
# using Pkg
# pkg "add OceanBioME, Oceananigans, CairoMakie, EnsembleKalmanProcesses, Distributions"
# ```

using OceanBioME, EnsembleKalmanProcesses, JLD2, CairoMakie, Oceananigans.Units, Oceananigans, BenchmarkTools

const year = years = 365day

# Setup the forward model

@inline PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

z = -10 # nominal depth of the box for the PAR profile
@inline PAR(t)::Float64 = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay 

function run_box_simulation()
    biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid = BoxModelGrid,
                                                                 light_attenuation_model = nothing)

    model = BoxModel(; biogeochemistry, forcing = (; PAR))

    set!(model, N = 10.0, P = 0.1, Z = 0.01)

    simulation = Simulation(model; Δt = 20minutes, stop_iteration = 1000, verbose = false)

    simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box_benchmarking.jld2", schedule = IterationInterval(20), overwrite_existing = true)

    @info "Running the model..."
    run!(simulation)
end

#####
##### results
#####

# Config: 1000 iterations with output every 8 hours and 20min steps
# origional @btime 317.5ms (1483607 allocations: 243.03 MiB)
# removing kernel launching from store tendencies: 291.546 ms (1293607 allocations: 187.95 MiB)
# removed kernel launching from rk3 substepping: 265.823 ms (1000607 allocations: 79.11 MiB)
# removed broadcasting in update state: 120.605 ms (619379 allocations: 63.98 MiB)
# no outputs:   23.523 ms (370344 allocations: 30.31 MiB)
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

using OceanBioME, EnsembleKalmanProcesses, JLD2, CairoMakie, Oceananigans.Units, Oceananigans
using LinearAlgebra, Random

using Distributions

using EnsembleKalmanProcesses
using EnsembleKalmanProcesses.ParameterDistributions

const year = years = 365day

rng_seed = 41
rng = Random.MersenneTwister(rng_seed)

# Setup the forward model

@inline PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

z = -10 # nominal depth of the box for the PAR profile
@inline PAR(t) = PAR⁰(t) * exp(0.2z) # Modify the PAR based on the nominal depth and exponential decay 

function run_box_simulation(initial_photosynthetic_slope,
                            base_maximum_growth,
                            nutrient_half_saturation,
                            phyto_base_mortality_rate,
                            j)

    biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid = BoxModelGrid,
                                                                 initial_photosynthetic_slope,
                                                                 base_maximum_growth,
                                                                 nutrient_half_saturation,
                                                                 phyto_base_mortality_rate,
                                                                 light_attenuation_model = nothing)

    model = BoxModel(; biogeochemistry, forcing = (; PAR))

    set!(model, N = 10.0, P = 0.1, Z = 0.01)

    simulation = Simulation(model; Δt = 20minutes, stop_time = 3years, verbose = false)

    simulation.output_writers[:fields] = JLD2OutputWriter(model, model.fields; filename = "box_calibration_$j.jld2", schedule = TimeInterval(8hours), overwrite_existing = true)

    @info "Running the model..."
    run!(simulation)

    P = FieldTimeSeries("box_calibration_$j.jld2", "P")

    times = P.times

    return P[1, 1, 1, length(times)-1092:end], times[length(times)-1092:end]
end

# Define the forward map

function G(u, j)
    (initial_photosynthetic_slope,
     base_maximum_growth,
     nutrient_half_saturation,
     phyto_base_mortality_rate) = u

    P, times = run_box_simulation(initial_photosynthetic_slope,
                                  base_maximum_growth,
                                  nutrient_half_saturation,
                                  phyto_base_mortality_rate,
                                  j)

    peak, winter, average, peak_timing, die_off_time = extract_observables(P, times)

    return [peak, winter, average, peak_timing, die_off_time], P
end

function extract_observables(P, times)
    if all(P .> 0) # model failure - including just in case
        peak = maximum(P)
        winter = minimum(P)
        average = mean(P)

        peak_timing = times[findmax(P)[2]]

        growth_rate = diff(P)[546:end]

        die_off_time = times[545 + findmin(growth_rate)[2]]

        return peak, winter, average, peak_timing./day, die_off_time./day
    else
        return NaN, NaN, NaN, NaN, NaN
    end
end

# Generate the "truth" data (normally you would load observations etc here)

Γ = Diagonal([0.001, 0.0001, 0.002, 5., 5.])

noise_dist = MvNormal(zeros(5), Γ)

truth = (0.15/day, 0.7/day, 2.4, 0.01/day)
obs, P₀ = G(truth, 1)

y = obs .+ rand(noise_dist)

# Solve the inverse problem and record all of the results for plotting purposes

prior_u1 = constrained_gaussian("initial_photosynthetic_slope", 0.1953 / day, 0.05 / day, 0, Inf)
prior_u2 = constrained_gaussian("base_maximum_growth", 0.6989 / day, 0.1/ day, 0, Inf)
prior_u3 = constrained_gaussian("nutrient_half_saturation", 2.3868, 0.5, 0, Inf)
prior_u4 = constrained_gaussian("phyto_base_mortality_rate", 0.0101 / day, 0.01 / day, 0, Inf)

prior = combine_distributions([prior_u1, prior_u2, prior_u3, prior_u4])

N_ensemble = 8
N_iterations = 5

initial_ensemble = construct_initial_ensemble(rng, prior, N_ensemble)

ensemble_kalman_process = EnsembleKalmanProcess(initial_ensemble, y, Γ, Inversion(); rng, failure_handler_method = SampleSuccGauss())

P = zeros(1093, N_ensemble, N_iterations) # recording all of the results for plotting only (not essential)

for i in 1:N_iterations
    @info "Iteration: $i"
    params_i = get_ϕ_final(prior, ensemble_kalman_process)

    G_ens = zeros(5, N_ensemble)

    Threads.@threads for j in 1:N_ensemble
        G_ens[:, j], P[:, j, i] = G(params_i[:, j], j)
    end

    update_ensemble!(ensemble_kalman_process, G_ens)
end

final_ensemble = get_ϕ_final(prior, ensemble_kalman_process)

# Plot the results

fig = Figure()

n = Observable(1)

title = @lift string("Generation: ", $n)

P_plts = [@lift P[:, j, $n] for j in 1:N_ensemble]

fig = Figure(size = (1200, 800));

ax = Axis(fig[1, 1], xlabel = "Day of year", ylabel = "Phytoplankton concentration (mmol/m³)"; title)

[lines!(ax, [1:8hours:365days-16hours;]./day, P_plts[j], color = :black, alpha = 0.2) for j in 1:N_ensemble]

lines!(ax, [1:8hours:365days-16hours;]./day, P₀, color = :black)

record(fig, "data_assimilation.mp4", 1:size(P, 3); framerate = 2) do i; n[] = i; end
nothing #hide

# ![](data_assimilation.mp4)

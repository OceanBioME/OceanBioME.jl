using OceanBioME, EnsembleKalmanProcesses, JLD2, CairoMakie, Oceananigans.Units
using LinearAlgebra, Random

using Distributions

using EnsembleKalmanProcesses
using EnsembleKalmanProcesses.ParameterDistributions

const year = years = 365day

rng_seed = 41
rng = Random.MersenneTwister(rng_seed)

#####
##### Setup model
#####

@inline PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline PAR(t) = PAR⁰(t) * exp(- 0.2 * 10)

function model(initial_photosynthetic_slope, 
               base_maximum_growth, 
               nutrient_half_saturation, 
               phyto_base_mortality_rate, 
               j)

    biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid = BoxModelGrid(),
                                                                 initial_photosynthetic_slope, 
                                                                 base_maximum_growth, 
                                                                 nutrient_half_saturation, 
                                                                 phyto_base_mortality_rate)

    model = BoxModel(; biogeochemistry, forcing = (; PAR))

    model.Δt = 20minutes
    model.stop_time = 5years

    set!(model, N = 10.0, P = 0.1, Z = 0.01)

    run!(model, save_interval = 24, feedback_interval = Inf, save = SaveBoxModel("box_calibration_$j.jld2"))

    file = jldopen("box_calibration_$j.jld2")
    times = parse.(Float64, keys(file["values"]))

    times = times[length(times)-1092:end]

    P = zeros(1093)

    # load a year of data
    for (idx, t_idx) in enumerate(length(times)-1092:length(times))
        P[idx] = file["values/$(times[t_idx])"].P
    end

    close(file)

    return P, times
end

function G(u, j)
    (initial_photosynthetic_slope, 
     base_maximum_growth, 
     nutrient_half_saturation, 
     phyto_base_mortality_rate) = u

    P, times = model(initial_photosynthetic_slope, 
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

#####
##### Generate the "truth" data (normally you would load observations etc here)
#####

# the "true" parameters
Γ = Diagonal([0.001, 0.0001, 0.002, 5., 5.])

noise_dist = MvNormal(zeros(5), Γ)

truth = (0.15/day, 0.7/day, 2.4, 0.01/day)
obs, P₀ = G(truth, 1)

y = obs  .+ rand(noise_dist)

#####
##### Solve the inverse problem
#####

# prior with default model values
prior_u1 = constrained_gaussian("initial_photosynthetic_slope", 0.1953 / day, 0.05 / day, 0, Inf)
prior_u2 = constrained_gaussian("base_maximum_growth", 0.6989 / day, 0.1/ day, 0, Inf)
prior_u3 = constrained_gaussian("nutrient_half_saturation", 2.3868, 0.5, 0, Inf)
prior_u4 = constrained_gaussian("phyto_base_mortality_rate", 0.0101 / day, 0.01 / day, 0, Inf)

prior = combine_distributions([prior_u1, prior_u2, prior_u3, prior_u4])

N_ensemble = 32
N_iterations = 20

initial_ensemble = construct_initial_ensemble(rng, prior, N_ensemble)

ensemble_kalman_process = EnsembleKalmanProcess(initial_ensemble, y, Γ, Inversion(); rng, failure_handler_method = SampleSuccGauss())

P = zeros(1093, N_ensemble, N_iterations) # saving all of the results for plotting purposes

for i in 1:N_iterations
    @info "Iteraition: $i"
    params_i = get_ϕ_final(prior, ensemble_kalman_process)

    G_ens = zeros(5, N_ensemble)

    Threads.@threads for j in 1:N_ensemble
        G_ens[:, j], P[:, j, i] = G(params_i[:, j], j)
    end

    update_ensemble!(ensemble_kalman_process, G_ens)
end

final_ensemble = get_ϕ_final(prior, ensemble_kalman_process)

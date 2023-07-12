dd = DataDep(
    "regression_results_npzd",
    "Base case for regression test", 
    "https://github.com/OceanBioME/OceanBioME_example_data/raw/main/regression_results/npzd.jld2"
)

register(dd)

function run_NPZD_box_regression_test(; comparison_path = datadep"regression_results_npzd/npzd.jld2")
    PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / 365days)) * (1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days)^2))) + 2

    z = -10
    PAR(t) = PAR⁰(t) * exp(0.2z)

    model = BoxModel(biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(grid = BoxModelGrid()), forcing = (; PAR))
    model.Δt = 5minutes
    model.stop_time = 2 * 10^3 * days

    set!(model, N = 10.0, P = 0.1, Z = 0.01)

    # ## Run the model (should only take a few seconds)
    @info "Running the model..."
    run!(model, save_interval = 100, save = SaveBoxModel("npzd.jld2"))

    @info "Model finished"

    vars = (:N, :P, :Z, :D, :PAR)

    times, timeseries = load_boxmodel("npzd.jld2", vars)
    _, comparison_timeseries = load_boxmodel(comparison_path, vars)

    @test all([all(timeseries[var] .≈ comparison_timeseries[var]) for var in vars])
end

run_NPZD_box_regression_test()
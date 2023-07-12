dd = DataDep(
    "regression_results_lobster_core",
    "Base case for regression test", 
    "https://github.com/OceanBioME/OceanBioME_example_data/raw/main/regression_results/lobster/core.jld2"
)

register(dd)

function run_LOBSTER_box_regression_test(; comparison_path = datadep"regression_results_lobster_core/core.jld2")
    PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / 365days)) * (1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days)^2))) + 2

    z = -10
    PAR(t) = PAR⁰(t) * exp(0.2z)

    model = BoxModel(biogeochemistry = LOBSTER(grid = BoxModelGrid()), forcing = (; PAR))
    model.Δt = 5minutes
    model.stop_time = 2 * 10^3 * days

    set!(model, NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01)

    # ## Run the model (should only take a few seconds)
    @info "Running the model..."
    run!(model, save_interval = 100, save = SaveBoxModel("lobster.jld2"))

    @info "Model finished"

    vars = (:NO₃, :NH₄, :P, :Z, :DOM, :sPOM, :bPOM, :PAR)

    times, timeseries = load_boxmodel("lobster.jld2", vars)
    _, comparison_timeseries = load_boxmodel(comparison_path, vars)

    @test all([all(timeseries[var] .≈ comparison_timeseries[var]) for var in vars])
end

dd = DataDep(
    "regression_results_lobster_carbon",
    "Base case for regression test", 
    "https://github.com/OceanBioME/OceanBioME_example_data/raw/main/regression_results/lobster/carbon.jld2"
)

register(dd)

function run_LOBSTER_carbon_box_regression_test(; comparison_path = datadep"regression_results_lobster_carbon/carbon.jld2")
    PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / 365days)) * (1 / (1 + 0.2 * exp(-((mod(t, 365days) - 200days) / 50days)^2))) + 2

    z = -10
    PAR(t) = PAR⁰(t) * exp(0.2z)

    model = BoxModel(biogeochemistry = LOBSTER(; grid = BoxModelGrid(), carbonates = true, oxygen = true, variable_redfield = true), forcing = (; PAR))
    model.Δt = 5minutes
    model.stop_time = 2 * 10^3 * days

    set!(model, NO₃ = 10.0, NH₄ = 0.1, P = 0.1, Z = 0.01, DIC = 2000, Alk = 1000, O₂ = 200)

    # ## Run the model (should only take a few seconds)
    @info "Running the model..."
    run!(model, save_interval = 100, save = SaveBoxModel("lobster_carbon.jld2"))

    @info "Model finished"

    vars = (:NO₃, :NH₄, :P, :Z, :DON, :sPON, :bPON, :DOC, :sPOC, :bPOC, :DIC, :Alk, :O₂)

    times, timeseries = load_boxmodel("lobster_carbon.jld2", vars)
    _, comparison_timeseries = load_boxmodel(comparison_path, vars)

    @test all([all(timeseries[var] .≈ comparison_timeseries[var]) for var in vars])
end

run_LOBSTER_box_regression_test()
run_LOBSTER_carbon_box_regression_test()
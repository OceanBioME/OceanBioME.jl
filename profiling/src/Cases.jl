module Cases
using Oceananigans
using OceanBioME
using Printf

using Oceananigans.Units
using Oceananigans.Fields: FunctionField, ConstantField

const year = years = 365days

"""
    simple_LOBSTER()

Run a simple representative case of 3D box with 🦞 LOBSTER 🦞 model

Based on the [single column example ](https://github.com/OceanBioME/OceanBioME.jl/blob/main/examples/column.jl)
slightly modified based on the [old benchmark case](https://github.com/OceanBioME/OceanBioME.jl/blob/main/benchmark/benchmark_simulations.jl#L12)
"""
function simple_LOBSTER(; backend = CPU(), fast_kill=false)

    # Parameters
    @inline PAR⁰(x, y, t) =
        60 *
        (1 - cos((t + 15days) * 2π / year)) *
        (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2
    @inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)
    @inline fmld1(t) =
        H(t, 50days, year) *
        (1 / (1 + exp(-(t - 100days) / 5days))) *
        (1 / (1 + exp((t - 330days) / 25days)))
    @inline MLD(t) = - (
        10 +
        340 * (
            1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))
        )
    )
    @inline κₜ(x, y, z, t) = 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4
    @inline temp(x, y, z, t) = 2.4 * cos(t * 2π / year + 50days) + 10
    grid_size = (4, 4, 50)

    # Setup case
    # TODO:
    #  - What about Time (T) and S fields from the example?
    grid = RectilinearGrid(;
        size = grid_size,
        extent = (20meters, 20meters, 200meters),
        topology = (Periodic, Periodic, Bounded),
    )

    biogeochemistry = LOBSTER(;
        grid,
        surface_photosynthetically_active_radiation = PAR⁰,
        carbonates = true,
        scale_negatives = true,
    )
    CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition()

    # Why do we need these?
    # TODO: Verify
    clock = Clock(; time = 0.0)
    T = FunctionField{Center,Center,Center}(temp, grid; clock)
    S = ConstantField(35.0)


    model = NonhydrostaticModel(;
        grid,
        closure = ScalarDiffusivity(ν = κₜ, κ = κₜ),
        biogeochemistry,
        boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),),
        auxiliary_fields = (; T, S),
    )
    set!(model, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2239.8, Alk = 2409.0)

    stop_time = fast_kill ? 3minutes : 10days;
    simulation = Simulation(model; Δt = 3minutes, stop_time)

    progress_message(sim) = @printf(
        "Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
        iteration(sim),
        prettytime(sim),
        prettytime(sim.Δt),
        prettytime(sim.run_wall_time)
    )
    add_callback!(simulation, progress_message, TimeInterval(1days))

    # Run
    run!(simulation)


    # Probably extract some results for regression check
    # TODO: Add this

end



end

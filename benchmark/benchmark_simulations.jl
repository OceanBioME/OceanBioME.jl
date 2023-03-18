using OceanBioME, Oceananigans, BenchmarkTools
using Oceananigans.Units
include("Benchmark.jl")

# Here I am attempting to replicate the setup of https://gmd.copernicus.org/articles/15/1567/2022/
# Although this is obviously very different (as the paper is about MPI) if our results are same
# order of mag it would be a good comparison to an established fast model
# I will compair to their Land Only Sub Domains Removed (LSR) cases but they do have bathymetry still
# at ~30 points there at about 1/2 hour per year of simulation - we're at about 10 mins on a single (vs ~4000) core so suspicious that this is not a valid comparison

function benchmark_config(n, biogeochemistry = nothing, callbacks = NamedTuple(), runtime = 5years)
    extent = n * 100kilometers
    grid = RectilinearGrid(size = (n, n, 16), extent = (extent, extent, 1kilometers)) 

    dTdz = 0.005 # K m⁻¹

    T_bcs = FieldBoundaryConditions(bottom = GradientBoundaryCondition(dTdz))

    u₁₀ = 0.1    # m s⁻¹, average wind velocity 10 meters above the ocean
    cᴰ = 2.5e-3 # dimensionless drag coefficient
    ρₐ = 1.225  # kg m⁻³, average density of air at sea-level
    ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean

    Qᵘ(x, y, t) = - ρₐ / ρₒ * cᴰ * u₁₀ * abs(u₁₀) * sin(π * x / 400kilometers)# m² s⁻²

    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))

    model = HydrostaticFreeSurfaceModel(; grid, biogeochemistry = biogeochemistry(; grid),
                                          coriolis = FPlane(f = 1e-4),
                                          boundary_conditions = (u = u_bcs, T = T_bcs),
                                          closure = ScalarDiffusivity(ν = 1e-5, κ = 1e-5))

    Ξ(z) = randn() * z/grid.Lz * (z/grid.Lz + 1)

    Tᵢ(x, y, z) = 20 + dTdz * z + dTdz * model.grid.Lz * 1e-6 * Ξ(z)

    uᵢ(x, y, z) = 1e-3 * Ξ(z)

    set!(model, u = uᵢ, v = uᵢ, S = 35, T = Tᵢ)

    Pᵢ(x, y, z) = rand()

    for (name, tracer) in pairs(model.tracers)
        if !(name in (:T, :S))
            set!(tracer, Pᵢ)
        end
    end

    # warm up
    time_step!(model, 1hour)

    benchmark = @benchmark begin
        time_step!($model, 1)
    end samples=10

    return benchmark
end

function benchmark_LOBSTER(n)
    grid = RectilinearGrid(size=(n, n, n), extent=(200, 200, 200)) 

    biogeochemistry = LOBSTER(; grid) 

    model = NonhydrostaticModel(; grid, biogeochemistry,
                                  advection = CenteredSecondOrder(),
                                  timestepper = :RungeKutta3,
                                  closure = ScalarDiffusivity(ν=1e-4, κ=1e-4))
                                  
    set!(model, P = 0.03, Z = 0.03, NO₃ = 11, NH₄ = 0.05)
    simulation = Simulation(model, Δt=10minutes, stop_iteration=1) 

    # warmup
    run!(simulation)

    trials = @benchmark run!($simulation) samples=5

    return trials
end


function benchmark_neg_protection(n)
    grid = RectilinearGrid(size=(n, n, n), extent=(200, 200, 200)) 

    biogeochemistry = LOBSTER(; grid) 

    model = NonhydrostaticModel(; grid, biogeochemistry,
                                  advection = CenteredSecondOrder(),
                                  timestepper = :RungeKutta3,
                                  closure = ScalarDiffusivity(ν=1e-4, κ=1e-4))
                                  
    set!(model, P = 0.03, Z = 0.03, NO₃ = 11, NH₄ = 0.05)
    simulation = Simulation(model, Δt=10minutes, stop_iteration=1) 

    scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM))
    simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())

    # warmup
    run!(simulation)

    trials = @benchmark run!($simulation) samples=5

    return trials
end

Ns = [2 ^ n for n in 2:6]

no_bgc(args...; kwargs...) = nothing

biogeochemistry = [no_bgc, LOBSTER]#, NutrientPhytoplanktonZooplanktonDetritus] # can't straightforwardly do this here because both NPZD model and buoyancy require T fields

suite = run_benchmarks(benchmark_config; Ns, biogeochemistry)

df = benchmarks_dataframe(suite)

benchmarks_pretty_table(df, title = "OceanBioME_benchmarks_8")
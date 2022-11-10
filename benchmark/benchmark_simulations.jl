using OceanBioME, Oceananigans, BenchmarkTools
using Oceananigans.Units
include("Benchmark.jl")

function benchmark_LOBSTER(n)
    grid = RectilinearGrid(size=(n, n, n), extent=(200, 200, 200)) 
    PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

    @inline PAR⁰(t) = 10.0

    bgc = Setup.Oceananigans(:LOBSTER, grid, LOBSTER.defaults) 

    model = NonhydrostaticModel(
        advection = WENO(;grid),
        timestepper = :RungeKutta3,
        grid = grid,
        tracers = bgc.tracers,
        closure = ScalarDiffusivity(ν=1e-4, κ=1e-4), 
        forcing = bgc.forcing,
        boundary_conditions = bgc.boundary_conditions,
        auxiliary_fields = (; PAR)
    )
    set!(model, P=0.03, Z=0.03, NO₃=11, NH₄=0.05)
    simulation = Simulation(model, Δt=10minutes, stop_iteration=1) 

    simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(LOBSTER.defaults, Light.twoBands.defaults), (surface_PAR=PAR⁰,)), TimeStepCallsite());
    simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))

    # warmup
    run!(simulation)

    trials = @benchmark run!($simulation) samples=5

    return trials
end

function benchmark_sediment(n)
    grid = RectilinearGrid(size=(n, n, n), extent=(200, 200, 200)) 
    PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

    @inline PAR⁰(t) = 10.0
    @inline t_function(x, y, z, t) = 13.0
    @inline s_function(x, y, z, t) = 35.0

    # Specify the boundary conditions for DIC and OXY based on the air-sea CO₂ and O₂ flux
    dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
    oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))

    # Have to prespecify fields for sediment model
    NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY = CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid), CenterField(grid)
    
    #sediment bcs
    sediment=Boundaries.Sediments.Soetaert.setupsediment(grid, (;D, DD); POM_w=(D=LOBSTER.D_sinking, DD=LOBSTER.DD_sinking))

    bgc = Setup.Oceananigans(:LOBSTER, grid, LOBSTER.defaults, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc), sinking=true, open_bottom=true, bottomboundaries=sediment.boundary_conditions)

    model = NonhydrostaticModel(
        advection = WENO(;grid),
        timestepper = :RungeKutta3,
        grid = grid,
        tracers = bgc.tracers,
        closure = ScalarDiffusivity(ν=1e-4, κ=1e-4), 
        forcing =  merge(bgc.forcing, sediment.forcing),
        boundary_conditions = bgc.boundary_conditions,
        auxiliary_fields = merge((; PAR), sediment.auxiliary_fields)
    )
    set!(model, P=0.03, Z=0.03, NO₃=11, NH₄=0.05)
    simulation = Simulation(model, Δt=10minutes, stop_iteration=1) 

    simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(LOBSTER.defaults, Light.twoBands.defaults), (surface_PAR=PAR⁰,)), TimeStepCallsite());
    simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))

    # warmup
    run!(simulation)

    trials = @benchmark run!($simulation) samples=5

    return trials
end

function benchmark_SLatissima(n)
    grid = RectilinearGrid(size=(n, n, n), extent=(200, 200, 200)) 
    PAR = Oceananigans.Fields.Field{Center, Center, Center}(grid)  

    @inline PAR⁰(t) = 10.0
    @inline t_function(x, y, z, t) = 13.0
    @inline s_function(x, y, z, t) = 35.0

    n_kelp = floor(Int, n/4) # number of kelp fronds
    z₀ = [-101+100/n_kelp:100/n_kelp:-1;]*1.0 # depth of kelp fronds

    kelp_particles = SLatissima.setup(n_kelp, Lx/2, Ly/2, z₀, 
                                                        30.0, 0.01, 0.1, 57.5;
                                                        scalefactor = 5.0, 
                                                        T = t_function, S = s_function, urel = 0.2, 
                                                        optional_tracers = (:NH₄, :DIC, :DD, :DDᶜ, :OXY, :DOM))

    dic_bc = Boundaries.airseasetup(:CO₂, forcings=(T=t_function, S=s_function))
    oxy_bc = Boundaries.airseasetup(:O₂, forcings=(T=t_function, S=s_function))
    bgc = Setup.Oceananigans(:LOBSTER, grid, LOBSTER.defaults, optional_sets=(:carbonates, :oxygen), topboundaries=(DIC=dic_bc, OXY=oxy_bc)) 

    model = NonhydrostaticModel(
        advection = WENO(;grid),
        timestepper = :RungeKutta3,
        grid = grid,
        tracers = bgc.tracers,
        closure = ScalarDiffusivity(ν=1e-4, κ=1e-4), 
        forcing = bgc.forcing,
        boundary_conditions = bgc.boundary_conditions,
        auxiliary_fields = (; PAR),
        particles = kelp_particles
    )
    set!(model, P=0.03, Z=0.03, NO₃=11, NH₄=0.05, DIC=2200.0, ALK=2400.0, OXY=240.0)
    simulation = Simulation(model, Δt=10minutes, stop_iteration=1) 

    simulation.callbacks[:update_par] = Callback(Light.twoBands.update!, IterationInterval(1), merge(merge(LOBSTER.defaults, Light.twoBands.defaults), (surface_PAR=PAR⁰,)), TimeStepCallsite());
    simulation.callbacks[:neg] = Callback(scale_negative_tracers!; parameters=(conserved_group=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM), warn=false))
    sim.callbacks[:couple_particles] = Callback(infinitesimal_particle_field_coupling!; callsite = TendencyCallsite())
    # warmup
    run!(simulation)

    trials = @benchmark run!($simulation) samples=5

    return trials
end

Ns = [16, 32]#, 64]

# Run and summarize benchmarks

for (model, name) in zip((:LOBSTER, :sediment, :SLatissima), ("LOBSTER", "sediment", "SLatissima"))
    @info "Benchmarking $name"
    benchmark_func = Symbol(:benchmark_, model)
    @eval begin
        suite = run_benchmarks($benchmark_func; Ns)
    end

    df = benchmarks_dataframe(suite)
    benchmarks_pretty_table(df, title=name * " benchmarks")
end
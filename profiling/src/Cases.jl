module Cases
using Oceananigans
using OceanBioME
using Random
using Printf


using Oceananigans.Units
using Oceananigans.Fields: FunctionField, ConstantField

import OceanBioME.Particles: fetch_output

const year = years = 365days

"""
    simple_LOBSTER()

Run a simple representative case of 3D box with 🦞 LOBSTER 🦞 model

Based on the [single column example ](https://github.com/OceanBioME/OceanBioME.jl/blob/main/examples/column.jl)
slightly modified based on the [old benchmark case](https://github.com/OceanBioME/OceanBioME.jl/blob/main/benchmark/benchmark_simulations.jl#L12)
"""
function simple_LOBSTER(; backend = CPU(), fast_kill = false)

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

"""
    big_LOBSTER()

A function version of the representative calculation script contributed in
04ec222a.

Stripped of any plotting and with option for fast exit for pre-compilation.

This script is intended as a profiling example for the LOBSTER biogeochemical model
and how it can be used with active particles to model the growth of sugar kelp.
It is intended as a typical use example for the OceanBioME.jl package.
"""
function big_LOBSTER(;
    backend = CPU(),
    grid_size :: Tuple{Int, Int, Int} = (32, 32, 8),
    n_particles :: Int = 5,
    fast_kill :: Bool = false,
    enable_io :: Bool = true,
    runlength_scale :: Float64 = 1.0,
    filename :: AbstractString = "LOBSTER",
    )

    duration = fast_kill ? 0.1day : 10days * runlength_scale # Duration of the simulation

    # set the number of gridpoints in each direction
    Nx, Ny, Nz = grid_size

    # set the domain size
    Lx = 1kilometer
    Ly = 1kilometer
    Lz = 100meters

    # Construct a grid with uniform grid spacing.
    grid = RectilinearGrid(backend, size = (Nx, Ny, Nz), extent = (Lx, Ly, Lz))

    # Set the Coriolis and buoyancy models.
    coriolis = FPlane(f = 1e-4) # [s⁻¹]
    buoyancy = SeawaterBuoyancy()

    # Specify parameters that are used to construct the background state.
    background_state_parameters = (
        M = 1e-4,       # s⁻¹, geostrophic shear
        f = coriolis.f, # s⁻¹, Coriolis parameter
        N = 1e-4,       # s⁻¹, buoyancy frequency
        H = grid.Lz,
        g = buoyancy.gravitational_acceleration,
        α = buoyancy.equation_of_state.thermal_expansion,
    )

    # We assume a background buoyancy ``B`` with a constant stratification and also a constant lateral
    # gradient (in the zonal direction). The background velocity components ``U`` and ``V`` are prescribed
    # so that the thermal wind relationship is satisfied, that is, ``f \partial_z U = - \partial_y B`` and
    # ``f \partial_z V = \partial_x B``.
    T(x, y, z, t, p) = (p.M^2 * x + p.N^2 * (z + p.H)) / (p.g * p.α)
    V(x, y, z, t, p) = p.M^2 / p.f * (z + p.H)

    V_field = BackgroundField(V, parameters = background_state_parameters)
    T_field = BackgroundField(T, parameters = background_state_parameters)

    # Specify some horizontal and vertical viscosity and diffusivity.
    νᵥ = κᵥ = 1e-4 # [m² s⁻¹]
    vertical_diffusivity = VerticalScalarDiffusivity(ν = νᵥ, κ = κᵥ)


    # ## Kelp Particle setup
    n = n_particles # number of kelp bundles
    z₀ =  LinRange(-21, -1, n) # depth of kelp fronds
    particles = SugarKelpParticles(n; grid, scalefactors = fill(2000/n, n)) # and we want them to look like there are 500 in each bundle
    # Initial conditions for the kelp particles.
    set!(particles, A = 10, N = 0.01, C = 0.1, z = z₀, x = Lx / 2, y = Ly / 2)

    # Setup the biogeochemical model with optional carbonate chemistry turned on.
    biogeochemistry = LOBSTER(;
        grid,
        carbonates = true,
        oxygen = true,
        variable_redfield = true,
        scale_negatives = true,
        open_bottom = true,
        particles,
    )

    DIC_bcs = FieldBoundaryConditions(top = CarbonDioxideGasExchangeBoundaryCondition())

    # Model instantiation
    model = NonhydrostaticModel(;
        grid,
        biogeochemistry,
        boundary_conditions = (DIC = DIC_bcs,),
        advection = WENO(),
        timestepper = :RungeKutta3,
        coriolis,
        tracers = (:T, :S),
        buoyancy,
        background_fields = (T = T_field, v = V_field),
        closure = vertical_diffusivity,
    )

    # ## Initial conditions
    # Start with a bit of random noise added to the background thermal wind and an arbitary
    # biogeochemical state.

    # Although not required, we also set the random seed to ensure reproducibility of the results.
    Random.seed!(11)

    Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

    Ũ = 1e-3
    uᵢ(x, y, z) = Ũ * Ξ(z)
    vᵢ(x, y, z) = Ũ * Ξ(z)

    set!(
        model,
        u = uᵢ,
        v = vᵢ,
        P = 0.03,
        Z = 0.03,
        NO₃ = 4.0,
        NH₄ = 0.05,
        DIC = 2200.0,
        Alk = 2409.0,
        S = 35,
        T = 20,
    )

    # ## Setup the simulation
    # Choose an appropriate initial timestep for this resolution and set up the simulation

    Δx = minimum_xspacing(grid, Center(), Center(), Center())
    Δy = minimum_yspacing(grid, Center(), Center(), Center())
    Δz = minimum_zspacing(grid, Center(), Center(), Center())

    Δt₀ = 0.75 * min(Δx, Δy, Δz) / V(0, 0, 0, 0, background_state_parameters)

    simulation = Simulation(model, Δt = Δt₀, stop_time = duration)

    # Adapt the time step while keeping the CFL number fixed.
    wizard = TimeStepWizard(cfl = 0.75, diffusive_cfl = 0.75, max_Δt = 30minutes)
    simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))
    nothing #hide

    # Create a progress message.
    progress(sim) = @printf(
        "i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
        sim.model.clock.iteration,
        prettytime(sim.model.clock.time),
        prettytime(sim.run_wall_time),
        prettytime(sim.Δt),
        AdvectiveCFL(sim.Δt)(sim.model)
    )

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))

    if enable_io
        # Here, we add some diagnostics to calculate and output.
        u, v, w = model.velocities # unpack velocity `Field`s

        # and also calculate the vertical vorticity.
        ζ = Field(∂x(v) - ∂y(u))

        # Periodically save the velocities and vorticity to a file.
        simulation.output_writers[:fields] = JLD2Writer(
            model,
            merge(model.tracers, (; u, v, w, ζ));
            schedule = TimeInterval(2hours),
            filename = "$(filename)_bgc.jld2",
            overwrite_existing = true,
        )

        function fetch_output(particles::BiogeochemicalParticles, model)
            return particles
        end

        simulation.output_writers[:particles] = JLD2Writer(
            model,
            (; model.biogeochemistry.particles),
            filename = "$(filename)_particles.jld2",
            schedule = TimeInterval(2hours),
            overwrite_existing = true,
        )
    end

    # Run the simulation
    run!(simulation)

end
end

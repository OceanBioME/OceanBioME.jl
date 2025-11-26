# This script is intended as a profiling example for the PISCES biogeochemical model 
# with MPI functionality

# ## Model setup
using OceanBioME, Oceananigans, Printf

using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units
using Random
using OceanBioME.Sediments: sinking_flux

using Plots # for plotting

const year = years = 365days
nothing #hide

# Set the architecture to use.
#arch = CPU()
# to use an NVIDIA GPU use the follwing lines:
using CUDA
arch = GPU()

plot_results = true # Set to `false` to not plot the results at the end of the simulation.

duration = 10days # Duration of the simulation

filename = "PISCES" # Base filename for output files

# set the number of gridpoints in each direction
Nx = 32
Ny = 32
Nz = 8

# set the domain size
Lx = 1kilometer
Ly = 1kilometer
Lz = 100meters

# Construct a grid with uniform grid spacing.
grid = RectilinearGrid(arch, size = (Nx, Ny, Nz), extent = (Lx, Ly, Lz))

clock = Clock(; time = 0.0)

# Set the Coriolis and buoyancy models.
coriolis = FPlane(f = 1e-4) # [s⁻¹]
buoyancy = SeawaterBuoyancy()

# Specify parameters that are used to construct the background state.
background_state_parameters = (M = 1e-4,       # s⁻¹, geostrophic shear
                               f = coriolis.f, # s⁻¹, Coriolis parameter
                               N = 1e-4,       # s⁻¹, buoyancy frequency
                               H = grid.Lz,
                               g = buoyancy.gravitational_acceleration,
                               α = buoyancy.equation_of_state.thermal_expansion)

# We assume a background buoyancy ``B`` with a constant stratification and also a constant lateral
# gradient (in the zonal direction). The background velocity components ``U`` and ``V`` are prescribed
# so that the thermal wind relationship is satisfied, that is, ``f \partial_z U = - \partial_y B`` and
# ``f \partial_z V = \partial_x B``.
T_func(x, y, z, t, p) = (p.M^2 * x + p.N^2 * (z + p.H)) / (p.g * p.α)
V_func(x, y, z, t, p) = p.M^2 / p.f * (z + p.H)

V_field = BackgroundField(V_func, parameters = background_state_parameters)
T_field = BackgroundField(T_func, parameters = background_state_parameters)

# Specify a vertical viscosity and diffusivity.
νᵥ = κᵥ = 1e-2 # [m² s⁻¹]
vertical_diffusivity = VerticalScalarDiffusivity(ν = νᵥ, κ = κᵥ)

# ## Surface PAR and turbulent vertical diffusivity based on idealised mixed layer depth
# Setting up idealised functions for PAR and diffusivity 

@inline PAR⁰(x, y, t) = 20.0

const MLD = -100.0

@inline κₜ(x, y, z, t) = (1e-2 * (1 + tanh((z - MLD) / 10.0)) / 2 + 1e-4)

κ_field = FunctionField{Center, Center, Center}(κₜ, grid; clock)

# Model setup
carbon_chemistry = CarbonChemistry()

biogeochemistry = PISCES(; grid, 
                           mixed_layer_depth = -100.0,
                           mean_mixed_layer_vertical_diffusivity = ConstantField(1e-2), # this is by default computed now
                           surface_photosynthetically_active_radiation = PAR⁰,
                           carbon_chemistry)

CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition(; carbon_chemistry)
O₂_flux = OxygenGasExchangeBoundaryCondition()

@info "Setting up the model..."
model = NonhydrostaticModel(; grid,
                                      advection = WENO(),
                                      buoyancy,
                                      coriolis,
                                      clock,
                                      closure = ScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ = κ_field),
                                      biogeochemistry,
                                      background_fields = (T = T_field, v = V_field))
                                      #boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux), O₂ = FieldBoundaryConditions(top = O₂_flux)))


@info "Setting initial values..."

# Although not required, we also set the random seed to ensure reproducibility of the results.
Random.seed!(11)

Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u = uᵢ, v=vᵢ, P = 0.1, PChl = 0.025, PFe = 0.005,
            D = 0.01, DChl = 0.003, DFe = 0.0006, DSi = 0.004,
            Z = 0.06, M = 0.5,
            DOC = 0.5,
            POC = 5.4, SFe = 0.34,
            GOC = 8.2, BFe = 0.5, PSi = 0.04, CaCO₃ = 10^-10,
            NO₃ = 10, NH₄ = 0.1, PO₄ = 5.0, Fe = 0.6, Si = 8.6,
            DIC = 2205, Alk = 2560, O₂ = 317, T=20, S = 35)

Δx = minimum_xspacing(grid, Center(), Center(), Center())
Δy = minimum_yspacing(grid, Center(), Center(), Center())
Δz = minimum_zspacing(grid, Center(), Center(), Center())

Δt₀ = 0.75 * min(Δx, Δy, Δz) / V_func(0, 0, 0, 0, background_state_parameters)

simulation = Simulation(model, Δt = Δt₀, stop_time = duration)

# Adapt the time step while keeping the CFL number fixed.
wizard = TimeStepWizard(cfl = 0.75, diffusive_cfl = 0.75, max_Δt = 30minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(5))

progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                iteration(sim),
                                prettytime(sim),
                                prettytime(sim.Δt),
                                prettytime(sim.run_wall_time))

add_callback!(simulation, progress_message, IterationInterval(1))

# Here, we add some diagnostics to calculate and output.
u, v, w = model.velocities # unpack velocity `Field`s

# and also calculate the vertical vorticity.
ζ = Field(∂x(v) - ∂y(u))

#add_callback!(simulation, update_temperature!, IterationInterval(1))
# Periodically save the tracers, velocities and vorticity to a file.
simulation.output_writers[:fields] = JLD2Writer(model, merge(model.tracers, (; u, v, w, ζ));
                                                      schedule = TimeInterval(1hour),
                                                      filename = "$(filename)_bgc.jld2",
                                                      overwrite_existing = true)
# ## Run!
# We are ready to run the simulation
run!(simulation)

if plot_results
    using JLD2, Plots, Statistics

    # Load the saved output,
    ζ = FieldTimeSeries("$(filename)_bgc.jld2", "ζ")
    P = FieldTimeSeries("$(filename)_bgc.jld2", "P")
    NO₃ = FieldTimeSeries("$(filename)_bgc.jld2", "NO₃")
    NH₄ = FieldTimeSeries("$(filename)_bgc.jld2", "NH₄")
    DIC = FieldTimeSeries("$(filename)_bgc.jld2", "DIC")
    Alk = FieldTimeSeries("$(filename)_bgc.jld2", "Alk")
    T = FieldTimeSeries("$(filename)_bgc.jld2", "T")
    S = FieldTimeSeries("$(filename)_bgc.jld2", "S")

    file = jldopen("$(filename)_bgc.jld2")
    iterations = keys(file["timeseries/t"])

    times = ζ.times

    xζ, yζ, zζ = nodes(ζ)
    xc, yc, zc = nodes(P)
    nothing #hide

    Pcolumn = zeros(length(zc), length(times))
    NO₃column = zeros(length(zc), length(times))
    DICcolumn = zeros(length(zc), length(times))
    NH₄column = zeros(length(zc), length(times))
    Alkcolumn = zeros(length(zc), length(times))
    Tcolumn = zeros(length(zc), length(times))
    Scolumn = zeros(length(zc), length(times))

    # loop over all saved times and make an animation
    anim = @animate for (i, iter) in enumerate(iterations[1:end])
        @printf("Plotting frame %d of %d\n", i, length(times))

        # extract a slice at the top of the domain
        ζₙ = interior(ζ[i], :, :, grid.Nz)
        Nₙ = interior(NO₃[i], :, :, grid.Nz) .+ interior(NH₄[i], :, :, grid.Nz)
        Pₙ = interior(  P[i], :, :, grid.Nz)
        DICₙ = interior(DIC[i], :, :, grid.Nz)

        # extract a column
        Pcolumn[:,i] = interior(P[i], 1, 1, :)
        NO₃column[:,i] = interior(NO₃[i], 1, 1, :)
        DICcolumn[:,i] = interior(DIC[i], 1, 1, :)
        NH₄column[:,i] = interior(NH₄[i], 1, 1, :)
        Alkcolumn[:,i] = interior(Alk[i], 1, 1, :)
        Tcolumn[:,i] = interior(T[i], 1, 1, :)
        Scolumn[:,i] = interior(S[i], 1, 1, :)

        lims = [(minimum(test), maximum(test)) for test in (  ζ[:, :, grid.Nz, :],
                                           NO₃[:, :, grid.Nz, :] .+ NH₄[:, :, grid.Nz, :],
                                             P[:, :, grid.Nz, :],
                                           DIC[:, :, grid.Nz, :])]

        kwargs = (xlabel = "x (m)", ylabel = "y (m)", aspect_ratio = :equal)

        hm1 = heatmap(xζ, yζ, ζₙ, colormap = :balance, colorrange = lims[1]; kwargs...)
        hm1_title = "Vorticity"

        hm2 = heatmap(xc, yc, Nₙ, colormap = :bamako, colorrange = lims[2]; kwargs...)
        hm2_title = "Nutrients"

        hm3 = heatmap(xc, yc, Pₙ, colormap = :batlow, colorrange = lims[3]; kwargs...)
        hm3_title = "Phytoplankton"

        hm4 = heatmap(xc, yc, DICₙ, colormap = :devon, colorrange = lims[4]; kwargs...)
        hm4_title = "DIC" 
        
        plot(hm1, hm2, hm3, hm4, layout = (2, 2), size = (1000, 800),
            title = [hm1_title, hm2_title, hm3_title, hm4_title])

end
nothing #hide

mp4(anim, "PISCES.mp4", fps = 15) # hide

end # end if plot_results


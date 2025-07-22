# This script is intended as a profiling example for the LOBSTER biogeochemical model
# and how it can be used with active particles to model the growth of sugar kelp.
# It is intended as a typical use example for the OceanBioME.jl package.

# ## Model setup
# We load the required packages. 
using OceanBioME, Oceananigans, Printf
using Oceananigans.Units
using Random
using JLD2
using Plots
using OceanBioME.Particles: BiogeochemicalParticles
import OceanBioME.Particles: required_particle_fields, required_tracers, coupled_tracers

# Set the architecture to use.
arch = CPU()
# to use an NVIDIA GPU use the follwing lines:
# using CUDA
# arch = GPU()

plot_results = true # Set to `false` to not plot the results at the end of the simulation.

duration = 10days # Duration of the simulation

filename = "LOBSTER" # Base filename for output files

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
T(x, y, z, t, p) = (p.M^2 * x + p.N^2 * (z + p.H)) / (p.g * p.α)
V(x, y, z, t, p) = p.M^2 / p.f * (z + p.H)

V_field = BackgroundField(V, parameters = background_state_parameters)
T_field = BackgroundField(T, parameters = background_state_parameters)

# Specify some horizontal and vertical viscosity and diffusivity.
νᵥ = κᵥ = 1e-4 # [m² s⁻¹]
vertical_diffusivity = VerticalScalarDiffusivity(ν = νᵥ, κ = κᵥ)


# ## Kelp Particle setup
n = 5 # number of kelp bundles
z₀ = [-21:5:-1;] * 1.0 # depth of kelp fronds
particles = SugarKelpParticles(n; grid, 
                               scalefactors = fill(2000, n)) # and we want them to look like there are 500 in each bundle
# Initial conditions for the kelp particles.
set!(particles, A = 10, N = 0.01, C = 0.1, z = z₀, x = Lx / 2, y = Ly / 2)

# Setup the biogeochemical model with optional carbonate chemistry turned on.
biogeochemistry = LOBSTER(; grid,
                            carbonates = true,
                            oxygen = true,
                            variable_redfield = true,
                            scale_negatives = true,
                            open_bottom = true,
                            particles)

DIC_bcs = FieldBoundaryConditions(top = CarbonDioxideGasExchangeBoundaryCondition())

# Model instantiation
model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              boundary_conditions = (DIC = DIC_bcs, ),
                              advection = WENO(grid),
                              timestepper = :RungeKutta3,
                              coriolis,
                              tracers = (:T, :S),
                              buoyancy,
                              background_fields = (T = T_field, v = V_field),
                              closure = vertical_diffusivity)

# ## Initial conditions
# Start with a bit of random noise added to the background thermal wind and an arbitary
# biogeochemical state.

# Although not required, we also set the random seed to ensure reproducibility of the results.
Random.seed!(11)

Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v=vᵢ, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2409.0, S = 35, T = 20)

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
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(sim.run_wall_time),
                        prettytime(sim.Δt),
                        AdvectiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))

# Here, we add some diagnostics to calculate and output.
u, v, w = model.velocities # unpack velocity `Field`s

# and also calculate the vertical vorticity.
ζ = Field(∂x(v) - ∂y(u))

# Periodically save the velocities and vorticity to a file.
simulation.output_writers[:fields] = JLD2Writer(model, merge(model.tracers, (; u, v, w, ζ));
                                                      schedule = TimeInterval(2hours),
                                                      filename = "$(filename)_bgc.jld2",
                                                      overwrite_existing = true)

import OceanBioME.Particles: fetch_output                                                      
function fetch_output(particles::BiogeochemicalParticles, model)
    return particles
end

simulation.output_writers[:particles] = JLD2Writer(model, (; model.biogeochemistry.particles),
                                                         filename = "$(filename)_particles.jld2",
                                                         schedule = TimeInterval(2hours),
                                                         overwrite_existing = true)                                                      

nothing #hide

# Run the simulation
run!(simulation)

if plot_results

    # Load the saved output,
    ζ = FieldTimeSeries("$(filename)_bgc.jld2", "ζ")
    P = FieldTimeSeries("$(filename)_bgc.jld2", "P")
    NO₃ = FieldTimeSeries("$(filename)_bgc.jld2", "NO₃")
    NH₄ = FieldTimeSeries("$(filename)_bgc.jld2", "NH₄")
    DIC = FieldTimeSeries("$(filename)_bgc.jld2", "DIC")

    file_particles = jldopen("$(filename)_particles.jld2")
    iterations = keys(file_particles["timeseries/t"])

    times = ζ.times

    xζ, yζ, zζ = nodes(ζ)
    xc, yc, zc = nodes(P)
    nothing #hide

    # loop over all saved times and make an animation
    anim = @animate for (i, iter) in enumerate(iterations[1:end])
        @printf("Plotting frame %d of %d\n", i, length(times))

        # extract a slice at the top of the domain
        ζₙ = interior(ζ[i], :, :, grid.Nz)
        Nₙ = interior(NO₃[i], :, :, grid.Nz) .+ interior(NH₄[i], :, :, grid.Nz)
        Pₙ = interior(  P[i], :, :, grid.Nz)
        DICₙ = interior(DIC[i], :, :, grid.Nz)

        # read in particle data at current timestep
        particles = file_particles["timeseries/particles/$iter"]

        lims = [(minimum(T), maximum(T)) for T in (  ζ[:, :, grid.Nz, :],
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

mp4(anim, "LOBSTER.mp4", fps = 15) # hide

end # end if plot_results

# ![](eady.mp4)

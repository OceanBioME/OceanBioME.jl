# # Biogeochemistry in submesoscale eddies in the Eady model
#
# In this example we setup a 3D model with a constant background buoyancy gradient with associated 
# thermal wind (the [Eady model](https://en.wikipedia.org/wiki/Eady_model)) with the [LOBSTER](@ref LOBSTER)
# biogeochemical model. This example demonstrates how to use biogeochemistry in a more complicated physical model.
# The parameters in this example correspond roughly to those used by [taylor2016](@citet) and result to the
# generation of a single submesoscale eddy.

# ## Install dependencies
# First we ensure we have the required dependencies installed
# ```julia
# using Pkg
# pkg "add OceanBioME, Oceananigans, CairoMakie"
# ```

# ## Model setup
# We load the required packages. Although not required, we also set the random seed to ensure
# reproducibility of the results.
using OceanBioME, Oceananigans, Printf
using Oceananigans.Units

using Random
Random.seed!(42)

# Construct a grid with uniform grid spacing.
grid = RectilinearGrid(size = (32, 32, 8), extent = (1kilometer, 1kilometer, 100meters))

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

# Setup the biogeochemical model with optional carbonate chemistry turned on.

biogeochemistry = LOBSTER(; grid,
                            carbonates = true,
                            open_bottom = true)

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

simulation = Simulation(model, Δt = Δt₀, stop_time = 10days)

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
simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.tracers, (; u, v, w, ζ));
                                                      schedule = TimeInterval(2hours),
                                                      filename = "eady_turbulence_bgc",
                                                      overwrite_existing = true)

nothing #hide

# Run the simulation
run!(simulation)

# Now load the saved output,
  ζ = FieldTimeSeries("eady_turbulence_bgc.jld2", "ζ")
  P = FieldTimeSeries("eady_turbulence_bgc.jld2", "P")
NO₃ = FieldTimeSeries("eady_turbulence_bgc.jld2", "NO₃")
NH₄ = FieldTimeSeries("eady_turbulence_bgc.jld2", "NH₄")
DIC = FieldTimeSeries("eady_turbulence_bgc.jld2", "DIC")

times = ζ.times

xζ, yζ, zζ = nodes(ζ)
xc, yc, zc = nodes(P)
nothing #hide

# and plot.

using CairoMakie

n = Observable(1)

  ζₙ = @lift interior(  ζ[$n], :, :, grid.Nz)
  Nₙ = @lift interior(NO₃[$n], :, :, grid.Nz) .+ interior(NH₄[$n], :, :, grid.Nz)
  Pₙ = @lift interior(  P[$n], :, :, grid.Nz)
DICₙ = @lift interior(DIC[$n], :, :, grid.Nz)

fig = Figure(size = (1600, 1600), fontsize = 20)

lims = [(minimum(T), maximum(T)) for T in (  ζ[:, :, grid.Nz, :],
                                           NO₃[:, :, grid.Nz, :] .+ NH₄[:, :, grid.Nz, :],
                                             P[:, :, grid.Nz, :],
                                           DIC[:, :, grid.Nz, :])]

axis_kwargs = (xlabel = "x (m)", ylabel = "y (m)", aspect = DataAspect())

ax1 = Axis(fig[1, 1]; title = "Vertical vorticity (1 / s)", axis_kwargs...)
hm1 = heatmap!(ax1, xζ, yζ, ζₙ, colormap = :balance, colorrange = lims[1])
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3]; title = "Nutrient (NO₃ + NH₄) concentration (mmol N / m³)", axis_kwargs...)
hm2 = heatmap!(ax2, xc, yc, Nₙ, colormap = Reverse(:bamako), colorrange = lims[2])
Colorbar(fig[1, 4], hm2)

ax3 = Axis(fig[2, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hm3 = heatmap!(ax3, xc, yc, Pₙ, colormap = Reverse(:batlow), colorrange = lims[3])
Colorbar(fig[2, 2], hm3)

ax4 = Axis(fig[2, 3]; title = "Dissolved inorganic carbon (mmol C / m³)", axis_kwargs...)
hm4 = heatmap!(ax4, xc, yc, DICₙ, colormap = Reverse(:devon), colorrange = lims[4])
Colorbar(fig[2, 4], hm4)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)

record(fig, "eady.mp4", 1:length(times), framerate = 12) do i
    n[] = i
end
nothing #hide

# ![](eady.mp4)

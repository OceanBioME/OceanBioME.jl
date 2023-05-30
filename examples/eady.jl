# # Biogeochemistry in submesoscale eddies in the Eady model
# In this example we will setup a 3D model with a constant background buoyancy gradient with associated 
# thermal wind (the [Eady model](https://en.wikipedia.org/wiki/Eady_model)) with the [LOBSTER](@ref LOBSTER) biogeochemical model. 
# This demonstrates how to use biogeochemistry in a more complicated physical model.
# The parameters correspond roughly to those in [taylor2016](@cite) and will generate a single submesoscale eddy

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg "add OceanBioME, Oceananigans, CairoMakie"
# ```

# ## Model setup
# First load the required packages
using OceanBioME, Oceananigans, Printf
using Oceananigans.Units

# Construct a grid with uniform grid spacing
Lz = 100meters
grid = RectilinearGrid(size=(32, 32, 8), extent = (1kilometer, 1kilometer, Lz))

# Set the Coriolis parameter
coriolis = FPlane(f = 1e-4) # [s⁻¹]

# Specify parameters that are used to construct the background state
background_state_parameters = ( M2 = 1e-8,       # s⁻¹, geostrophic shear
                                 f = coriolis.f, # s⁻¹, Coriolis parameter
                                 N = 1e-4,       # s⁻¹, buoyancy frequency
                                Lz = Lz )

# Here, ``B`` is the background buoyancy field and ``V`` is the corresponding thermal wind
V(x, y, z, t, p) = p.M2 / p.f * (z - p.Lz/2)
B(x, y, z, t, p) = p.M2 * x + p.N^2 * (z - p.Lz/2)

V_field = BackgroundField(V, parameters = background_state_parameters)
B_field = BackgroundField(B, parameters = background_state_parameters)

# Specify some horizontal and vertical viscosity/diffusivity
νₕ = κₕ = 1e2  # [m² s⁻¹]
horizontal_diffusivity = HorizontalScalarDiffusivity(ν = νₕ, κ = κₕ)
#-
νᵥ = κᵥ = 1e-4 # [m² s⁻¹]
vertical_diffusivity = VerticalScalarDiffusivity(ν = νᵥ, κ = κᵥ)

# Setup the biogeochemical model with optional carbonate chemistry turned on

biogeochemistry = LOBSTER(; grid, 
                            carbonates = true,
                            open_bottom = true)

# To-do: change to a buoyancy parameterisation so we don't have to fake the temperature and salinity
DIC_bcs = FieldBoundaryConditions(top = GasExchange(; gas = :CO₂, temperature = (args...) -> 12, salinity = (args...) -> 35))

# Model instantiation
model = NonhydrostaticModel(; grid,
                              biogeochemistry,
                              boundary_conditions = (DIC = DIC_bcs, ),
                              advection = WENO(grid),
                              timestepper = :RungeKutta3,
                              coriolis,
                              tracers = :b,
                              buoyancy = BuoyancyTracer(),
                              background_fields = (b = B_field, v = V_field),
                              closure = (horizontal_diffusivity, vertical_diffusivity))

# ## Initial conditions
# Start with a bit of random noise added to the background thermal wind and an arbitary biogeochemical state
Ξ(z) = randn() * z / grid.Lz * (z / grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v=vᵢ, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2409.0)

simulation = Simulation(model, Δt = 15minutes, stop_time = 10days)

# Adapt the time step while keeping the CFL number fixed
wizard = TimeStepWizard(cfl=0.6, max_change = 1.5, max_Δt = 30minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(1))
nothing #hide

# Create a progress message 
progress(sim) = @printf("i: % 6d, sim time: % 10s, wall time: % 10s, Δt: % 10s, CFL: %.2e\n",
                        sim.model.clock.iteration,
                        prettytime(sim.model.clock.time),
                        prettytime(sim.run_wall_time),
                        prettytime(sim.Δt),
                        AdvectiveCFL(sim.Δt)(sim.model))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(10))

# Here, add some diagnostics to calculate and output

u, v, w = model.velocities # unpack velocity `Field`s

# and also calculate the vertical vorticity [s⁻¹].
ζ = Field(∂x(v) - ∂y(u))

# Periodically write the velocity and vorticity out to a file
simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.tracers, (; u, v, w, ζ));
                                                      schedule = TimeInterval(2hours),
                                                      filename = "eady_turbulence_bgc",
                                                      overwrite_existing = true)

# Prevent the tracer values going negative - this is especially important in this model while no positivity preserving diffusion is implimented
scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())
simulation.callbacks[:nan_tendencies] = Callback(remove_NaN_tendencies!; callsite = TendencyCallsite())
simulation.callbacks[:abort_zeros] = Callback(zero_negative_tracers!; callsite = UpdateStateCallsite())
nothing #hide

# Run the simulation
run!(simulation)

# Now load the saved output
  ζ = FieldTimeSeries("eady_turbulence_bgc.jld2", "ζ")
  P = FieldTimeSeries("eady_turbulence_bgc.jld2", "P")
NO₃ = FieldTimeSeries("eady_turbulence_bgc.jld2", "NO₃")
NH₄ = FieldTimeSeries("eady_turbulence_bgc.jld2", "NH₄")
DIC = FieldTimeSeries("eady_turbulence_bgc.jld2", "DIC")

times = ζ.times

xζ, yζ, zζ = nodes(ζ)
xc, yc, zc = nodes(P)

# and plot

using CairoMakie

n = Observable(1)

  ζₙ = @lift interior(  ζ[$n], :, :, grid.Nz)
  Nₙ = @lift interior(NO₃[$n], :, :, grid.Nz) .+ interior(NH₄[$n], :, :, grid.Nz)
  Pₙ = @lift interior(  P[$n], :, :, grid.Nz)
DICₙ = @lift interior(DIC[$n], :, :, grid.Nz)

fig = Figure(resolution = (1600, 1600), fontsize = 20)

lims = [(minimum(T), maximum(T)) for T in (  ζ[:, :, grid.Nz, :],
                                           NO₃[:, :, grid.Nz, :] .+ NH₄[:, :, grid.Nz, :],
                                             P[:, :, grid.Nz, :],
                                           DIC[:, :, grid.Nz, :])]

axis_kwargs = (xlabel = "x (m)", ylabel = "y (m)", aspect = DataAspect())

ax1 = Axis(fig[1, 1]; title = "Vertical vorticity (1 / s)", axis_kwargs...)
hm1 = heatmap!(ax1, xζ, yζ, ζₙ, levels = 33, colormap = :balance, colorrange = lims[1])
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3]; title = "Nutrient (NO₃ + NH₄) concentration (mmol N / m³)", axis_kwargs...)
hm2 = heatmap!(ax2, xc, yc, Nₙ, levels = 33, colormap = Reverse(:bamako), colorrange = lims[2])
Colorbar(fig[1, 4], hm2)

ax3 = Axis(fig[2, 1]; title = "Phytoplankton concentration (mmol N / m³)", axis_kwargs...)
hm3 = heatmap!(ax3, xc, yc, Pₙ, levels = 33, colormap = Reverse(:batlow), colorrange = lims[3])
Colorbar(fig[2, 2], hm3)

ax4 = Axis(fig[2, 3]; title = "Dissolved inorganic carbon (mmol C / m³)", axis_kwargs...)
hm4 = heatmap!(ax4, xc, yc, DICₙ, levels = 33, colormap = Reverse(:devon), colorrange = lims[4])
Colorbar(fig[2, 4], hm4)

title = @lift "t = $(prettytime(times[$n]))"
Label(fig[0, :], title, fontsize = 30)

record(fig, "eady.mp4", 1:length(times), framerate = 12) do i
    @info string("Plotting frame ", i, " of ", length(times))
    n[] = i
end
nothing #hide

# ![](eady.mp4)

# # Biogeochemistry in submesoscale eddies in the Eady model
# In this example we will setup a 3D model with a constant background buoyancy gradient with associated 
# thermal wind (the [Eady model](https://en.wikipedia.org/wiki/Eady_model)) with the [LOBSTER](@ref LOBSTER) biogeochemical model. 
# This demonstrates how to use biogeochemistry in a more complicated physical model.
# The parameters correspond roughly to those in [taylor2016](@cite) and will generate a single submesoscale eddy

# ## Install dependencies
# First we will check we have the dependencies installed
# ```julia
# using Pkg
# pkg "add OceanBioME, Oceananigans, Printf, CairoMakie"
# ```

# ## Model setup
# First load the required packages
using OceanBioME, Oceananigans, Printf
using Oceananigans.Units

# Construct a grid with uniform grid spacing
Lz = 100
grid = RectilinearGrid(size=(32, 32, 8), extent = (1kilometer, 1kilometer, Lz))

# Set the Coriolis parameter
coriolis = FPlane(f=1e-4) # [s⁻¹]

# Specify parameters that are used to construct the background state
background_state_parameters = ( M2 = 1e-8, # s⁻¹, geostrophic shear
                                f = coriolis.f,      # s⁻¹, Coriolis parameter
                                N = 1e-4)            # s⁻¹, buoyancy frequency

# Here, B is the background buoyancy field and V is the corresponding thermal wind
V(x, y, z, t, p) = + p.M2/p.f * (z - Lz/2)
B(x, y, z, t, p) = p.M2 * x + p.N^2 * (z - Lz/2)

V_field = BackgroundField(V, parameters = background_state_parameters)
B_field = BackgroundField(B, parameters = background_state_parameters)

# Specify the horizontal and vertical viscosity/diffusivity
κ₂z = 1e-4 # [m² s⁻¹] Vertical vertical viscosity and diffusivity
κ₂h = 1e-2 # [m² s⁻¹] Horizontal viscosity and diffusivity

vertical_diffusivity = VerticalScalarDiffusivity(ν=κ₂z, κ=κ₂z)
horizontal_diffusivity = HorizontalScalarDiffusivity(ν=κ₂h, κ=κ₂h)

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
                              advection = CenteredSecondOrder(),
                              timestepper = :RungeKutta3,
                              coriolis,
                              tracers = :b,
                              buoyancy = BuoyancyTracer(),
                              background_fields = (b = B_field, v = V_field),
                              closure = (vertical_diffusivity, horizontal_diffusivity))

# ## Initial conditions
# Start with a bit of random noise added to the background thermal wind and an arbitary biogeochemical state
Ξ(z) = randn() * z/grid.Lz * (z/grid.Lz + 1)

Ũ = 1e-3
uᵢ(x, y, z) = Ũ * Ξ(z)
vᵢ(x, y, z) = Ũ * Ξ(z)

set!(model, u=uᵢ, v=vᵢ, P = 0.03, Z = 0.03, NO₃ = 4.0, NH₄ = 0.05, DIC = 2200.0, Alk = 2409.0)

simulation = Simulation(model, Δt = 15minutes, stop_time = 10days)

# Adapt the time step while keeping the CFL number fixed
wizard = TimeStepWizard(cfl=0.85, max_change = 1.5, max_Δt = 10minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

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

# calculate the vertical vorticity [s⁻¹]
ζ = Field(∂x(v) - ∂y(u))

# horizontal divergence [s⁻¹]
δ = Field(∂x(u) + ∂y(v))

# Periodically write the velocity, vorticity, and divergence out to a file
simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.tracers, (; u, v, w, ζ, δ));
                                                      schedule = TimeInterval(4hours),
                                                      filename = "eady_turbulence_bgc",
                                                      overwrite_existing = true);

# Prevent the tracer values going negative - this is especially important in this model while no positivity preserving diffusion is implimented
scale_negative_tracers = ScaleNegativeTracers(; model, tracers = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM))
simulation.callbacks[:neg] = Callback(scale_negative_tracers; callsite = UpdateStateCallsite())
simulation.callbacks[:nan_tendencies] = Callback(remove_NaN_tendencies!; callsite = TendencyCallsite())
simulation.callbacks[:abort_zeros] = Callback(zero_negative_tracers!; callsite = UpdateStateCallsite())

# Run the simulation
run!(simulation)

# Now plot the results
using CairoMakie, JLD2

# Open the file with our data
file = jldopen(simulation.output_writers[:fields].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

times = zeros(length(iterations))

w, P, N, DIC = ntuple(n -> zeros(grid.Nx, grid.Ny, grid.Nz, length(iterations)), 5);

for (idx, it) in enumerate(iterations)
  w[:, :, :, idx] = file["timeseries/w/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz]
  P[:, :, :, idx] = file["timeseries/P/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz]
  N[:, :, :, idx] = file["timeseries/NO₃/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz] .+ file["timeseries/NH₄/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz]
  DIC[:, :, :, idx] = file["timeseries/DIC/$it"][1:grid.Nx, 1:grid.Ny, 1:grid.Nz]

  times[idx] = file["timeseries/t/$it"]
end

# Plot

fig = Figure(resolution = (1600, 1600))

n = Observable(1)

w_plt = @lift w[:, :, grid.Nz, $n]
N_plt = @lift N[:, :, grid.Nz, $n]
P_plt = @lift P[:, :, grid.Nz, $n]
DIC_plt = @lift DIC[:, :, grid.Nz, $n]

xs = xnodes(grid, Center(), Center(), Center())[1:grid.Nx]
ys = ynodes(grid, Center(), Center(), Center())[1:grid.Ny]
zs = znodes(grid, Center(), Center(), Center())[1:grid.Nz]

lims = [(minimum(T), maximum(T)) for T in (w, N, P, DIC)]

ax1 = Axis(fig[1, 1], aspect = DataAspect(), title = "Vertical velocity (m / s)")
hm1 = heatmap!(ax1, xs, ys, w_plt, levels = 33, colormap = :lajolla, colorrange = lims[1], interpolate = true)
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3], aspect = DataAspect(), title = "Nutrient (NO₃ + NH₄) concentration (mmol N / m³)")
hm2 = heatmap!(ax2, xs, ys, N_plt, levels = 33, colormap = Reverse(:bamako), colorrange = lims[2], interpolate = true)
Colorbar(fig[1, 4], hm2)

ax3 = Axis(fig[2, 1], aspect = DataAspect(), title = "Phytoplankton concentration (mmol N / m³)")
hm3 = heatmap!(ax3, xs, ys, P_plt, levels = 33, colormap = Reverse(:batlow), colorrange = lims[3], interpolate = true)
Colorbar(fig[2, 2], hm3)

ax4 = Axis(fig[2, 3], aspect = DataAspect(), title = "Dissolved inorganic carbon (mmol C / m³)")
hm4 = heatmap!(ax4, xs, ys, DIC_plt, levels = 33, colormap = Reverse(:devon), colorrange = lims[4], interpolate = true)
Colorbar(fig[2, 4], hm4)

supertitle = Label(fig[0, :], "t = 0.0", fontsize = 30)

record(fig, "eady.gif", 1:length(times), framerate = 10) do i
    n[] = i
    msg = string("Plotting frame ", i, " of ", length(times))
    print(msg * " \r")
    supertitle.text = "t=$(prettytime(times[i]))"
end

# ![Results](eady.gif)

using Oceananigans
using Oceananigans.Units
using OceanBioME.NPZD: NutrientPhytoplanktonZooplanktonDetritus
using OceanBioME.Light: twoBands
using OceanBioME: LOBSTER

grid = RectilinearGrid(size=(32, 32, 32), extent=(2000, 2000, 64), topology=(Periodic, Periodic, Bounded))
PAR = Oceananigans.Fields.CenterField(grid)

Qʰ = 200.0  # W m⁻², surface _heat_ flux
ρₒ = 1026.0 # kg m⁻³, average density at the surface of the world ocean
cᴾ = 3991.0 # J K⁻¹ kg⁻¹, typical heat capacity for seawater

Qᵀ(x, y, t) = Qʰ* ifelse(t>5days, -max(0.0, sin(t*π/12hours)), 0.0) / (ρₒ * cᴾ) # K m s⁻¹, surface _temperature_ flux
dTdz = 0.01 # K m⁻¹

T_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵀ),
                                bottom = GradientBoundaryCondition(dTdz))

u₁₀ = 10    # m s⁻¹, average wind velocity 10 meters above the ocean
cᴰ = 2.5e-3 # dimensionless drag coefficient
ρₐ = 1.225  # kg m⁻³, average density of air at sea-level

Qᵘ(x, y, t) = ifelse(t>5days, - ρₐ / ρₒ * cᴰ * u₁ ₀ * abs(u₁₀), 0.0) # m² s⁻²
u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵘ))
@inline Qˢ(x, y, t, S, evaporation_rate) = - evaporation_rate(t) * S # [salinity unit] m s⁻¹
evaporation_rate(t) = ifelse(t>5days, 1e-3 / hour * max(0.0, sin(t*π/12hours)) + 1e-4/hour, 0.0)# m s⁻¹
evaporation_bc = FluxBoundaryCondition(Qˢ, field_dependencies=:S, parameters=evaporation_rate)
S_bcs = FieldBoundaryConditions(top=evaporation_bc)
model = NonhydrostaticModel(; grid,
                            advection = UpwindBiasedFifthOrder(),
                            timestepper = :RungeKutta3,
                            closure = ScalarDiffusivity(ν=1e-4, κ=1e-4),
                            coriolis = FPlane(f=1e-4),
                            tracers = (:S, :T, :N, :P, :Z, :D), # P for Plankton
                            buoyancy = SeawaterBuoyancy(),
                            boundary_conditions = (u=u_bcs, T=T_bcs, S=S_bcs),
                            biogeochemistry=NutrientPhytoplanktonZooplanktonDetritus(),
                            auxiliary_fields = (; PAR))

set!(model, S=35.0, u=1.0, T=15.0, N=10.0, P=0.01, Z=0.01, D=0.0)


simulation = Simulation(model, Δt=0.5minutes, stop_time=10day)
Rd_chl = LOBSTER.defaults.Rd_chl
surface_PAR(t) = ifelse(t > 5days, 50 * max(0.0, sin(t*π/12hours)), 25)
simulation.callbacks[:update_par] = Callback(twoBands.update!, IterationInterval(1), merge(twoBands.defaults, (;surface_PAR, Rd_chl)), TimeStepCallsite());

wizard = TimeStepWizard(cfl=0.5, diffusive_cfl=0.5, max_Δt = 1minutes)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
using Printf

progress(sim) = @printf("Iteration: %d, time: %s, Δt: %s\n",
                        iteration(sim), prettytime(sim), prettytime(sim.Δt))

simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

# and a basic `JLD2OutputWriter` that writes velocities and both
# the two-dimensional and horizontally-averaged plankton concentration,

simulation.output_writers[:simple_output] =
    JLD2OutputWriter(model, merge(model.velocities, model.tracers),
                     schedule = TimeInterval(20minutes),
                     filename = "convecting_plankton.jld2",
                     overwrite_existing = true)

# !!! info "Using multiple output writers"
#     Because each output writer is associated with a single output `schedule`,
#     it often makes sense to use _different_ output writers for different types of output.
#     For example, smaller outputs that consume less disk space may be written more
#     frequently without threatening the capacity of your hard drive.
#     An arbitrary number of output writers may be added to `simulation.output_writers`.
#
# The simulation is set up. Let there be plankton:

run!(simulation)

# Notice how the time-step is reduced at early times, when turbulence is strong,
# and increases again towards the end of the simulation when turbulence fades.

# ## Visualizing the solution
#
# We'd like to a make a plankton movie. First we load the output file
# and build a time-series of the buoyancy flux,

filepath = simulation.output_writers[:simple_output].filepath

w_timeseries = FieldTimeSeries(filepath, "u")
T_timeseries = FieldTimeSeries(filepath, "T")
N_timeseries = FieldTimeSeries(filepath, "N")
P_timeseries = FieldTimeSeries(filepath, "P")
Z_timeseries = FieldTimeSeries(filepath, "Z")
D_timeseries = FieldTimeSeries(filepath, "D")

times = w_timeseries.times
# and then we construct the ``x, z`` grid,

w_lim = maximum(abs, interior(w_timeseries))
w_lims = (-w_lim, w_lim)
T_lims = (minimum(T_timeseries), maximum(T_timeseries))
N_lims = (minimum(N_timeseries), maximum(N_timeseries))
P_lims = (minimum(P_timeseries), maximum(P_timeseries))
Z_lims = (minimum(Z_timeseries), maximum(Z_timeseries))
D_lims = (minimum(D_timeseries), maximum(D_timeseries))

xw, yw, zw = nodes(w_timeseries)
xT, yT, zT = nodes(T_timeseries)

# Finally, we animate plankton mixing and blooming,

using GLMakie

@info "Making a movie about plankton..."

n = Observable(1)

title = @lift @sprintf("t = %s", prettytime(times[$n]))

wₙ = @lift interior(w_timeseries[$n], :, :, grid.Nz)
Tₙ = @lift interior(T_timeseries[$n], :, 1, :)
Nₙ = @lift interior(N_timeseries[$n], :, 1, :)
Pₙ = @lift interior(P_timeseries[$n], :, 1, :)
Zₙ = @lift interior(Z_timeseries[$n], :, 1, :)
Dₙ = @lift interior(D_timeseries[$n], :, 1, :)

fig = Figure(resolution = (1200, 2000))

ax_w = Axis(fig[2, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_T = Axis(fig[3, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_N = Axis(fig[4, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_P = Axis(fig[5, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_Z = Axis(fig[6, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)
ax_D = Axis(fig[7, 2]; xlabel = "x (m)", ylabel = "z (m)", aspect = 1)

fig[1, 1:2] = Label(fig, title, tellwidth=false)

hm_w = heatmap!(ax_w, xw, yw, wₙ; colormap = :balance, colorrange = w_lims)
Colorbar(fig[2, 1], hm_w; label = "x velocity (m s⁻¹)", flipaxis = false)
hm_T = heatmap!(ax_T, xT, zT, Tₙ; colormap = :balance, colorrange = T_lims)
Colorbar(fig[3, 1], hm_T; label = "Temperature (°C)", flipaxis = false)
hm_N = heatmap!(ax_N, xw, zw, Nₙ; colormap = :balance, colorrange = N_lims)
Colorbar(fig[4, 1], hm_N; label = "Nutrient Concentration (mmol N/m³)", flipaxis = false)
hm_P = heatmap!(ax_P, xw, zw, Pₙ; colormap = :balance, colorrange = P_lims)
Colorbar(fig[5, 1], hm_P; label = "Phytoplankton Concentration (mmol N/m³)", flipaxis = false)
hm_Z = heatmap!(ax_Z, xw, zw, Zₙ; colormap = :balance, colorrange = Z_lims)
Colorbar(fig[6, 1], hm_Z; label = "Zooplankton Concentration (mmol N/m³)", flipaxis = false)
hm_D = heatmap!(ax_D, xw, zw, Dₙ; colormap = :balance, colorrange = D_lims)
Colorbar(fig[7, 1], hm_D; label = "Detritus Concentration (mmol N/m³)", flipaxis = false)
# And, finally, we record a movie.

frames = 1:length(times)

@info "Making an animation of convecting plankton..."

record(fig, "convecting_plankton.mp4", frames, framerate=8) do i
    msg = string("Plotting frame ", i, " of ", frames[end])
    print(msg * " \r")
    n[] = i
end
nothing #hide

# ![](convecting_plankton.mp4)

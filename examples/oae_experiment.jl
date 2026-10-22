# # Idealised ocean alkalinity enhancement (OAE)
#
# In this example we setup a 3D model driven by a constant wind stress, with 3 replicates
# of the carbon chemistry system (DIC, alkalinity, and particulate inorganic carbon), labelled:
#
# 1. `control`, which receives no alkalinity,
# 2. `release`, which receives the alkalinity release, with secondary precipitation active,
# 3. `no_precipitation`, which receives the same release with secondary precipitation disabled
#    (`calcium_carbonate_precipitation_rate = 0`).
#
# The `release` and `control` replicates differ only by the alkalinity release, so their difference
# is the OAE signal; the `release` and `no_precipitation` replicates differ only by the precipitation
# rate, so their difference is what secondary precipitation costs.
#
# We use [`ExplicitCalciumCarbonate`](@ref) rather than [`CarbonateSystem`](@ref) so that the
# particulate inorganic carbon formed by the calcifiers, and by abiotic precipitation once the
# alkalinity addition pushes the saturation state `Ω` up, is carried as its own sinking `CaCO₃`
# tracer rather than being routed straight back into `DIC`/`Alk`.

# ## Install dependencies
# First we ensure we have the required dependencies installed
# ```julia
# using Pkg
# pkg "add OceanBioME, Oceananigans, CairoMakie"
# ```

# ## Model setup

using Oceananigans, OceanBioME, Oceananigans.Units, CUDA

using OceanBioME.Models.CarbonChemistryModel: KSP_aragonite

# we define the grid and the biogeochemistry, adding in the `inorganic_carbon` with 3 labelled replicates
grid = RectilinearGrid(GPU();
                       size = (64, 64, 8),
                       extent = (500, 500, 15))

precipitation_rate = 0.03 / day

# use aragonite rather than the default calcite solubility, as that is the phase usually
# discussed in the OAE secondary precipitation literature
carbon_chemistry = CarbonChemistry(; calcium_carbonate_solubility = KSP_aragonite())

biogeochemistry = LOBSTER(grid;
                          inorganic_carbon =
                              ExplicitCalciumCarbonate(grid;
                                                       replicates = (:control, :release, :no_precipitation),
                                                       carbon_chemistry,
                                                       calcium_carbonate_precipitation_rate =
                                                           (control = precipitation_rate,
                                                            release = precipitation_rate,
                                                            no_precipitation = 0),
                                                       calcium_carbonate_precipitation_exponent = 3))

# Secondary precipitation is off by default, so we turn it on for the `control` and `release`
# replicates. The rate law is
#
# ```math
# \text{precipitation} = k \max(0, \Omega - 1)^n [CaCO_3],
# ```
#
# following the empirical description of calcite precipitation from seawater of [Zuddas1994](@citet)
# and [Zuddas1998](@citet). This example is a fictional scenario with fictional parameter values for
# `k` and `n`, chosen only to give a visible, numerically well-behaved precipitation signal - they
# should not be read as literature values or as a prediction of a real deployment.

# Each replicate's tracers are named after its label (e.g. `DIC_control`, `Alk_control`, `CaCO₃_control`,
# and the saturation state `Ω_control`), which `carbonate_replicates` lists
replicates = carbonate_replicates(biogeochemistry)

# Next we construct the boundary conditions: the air-sea CO₂ exchange of each replicate, computed
# from its own DIC and alkalinity, and the wind stress
wind_speed = 5 # m/s
Cᴰ = 2e-3
ρₐ = 1.2
ρₒ = 1026

wind = FluxBoundaryCondition(-ρₐ/ρₒ * Cᴰ * wind_speed^2)

boundary_conditions = merge(carbon_dioxide_gas_exchange_boundary_conditions(biogeochemistry; wind_speed),
                            (; u = FieldBoundaryConditions(top = wind)))

# Next we define an alkalinity release in a circle in the center for 1 hour

@inline alkalinity_release(x, y, z, t, params) =
    ifelse((params.start_time <= t < params.start_time + params.duration) & 
            (z > params.depth) & 
            ((x-250)^2 + (y-250)^2 < params.radius^2), 
            params.release_rate, 0)

total_release = 2.5e7 # 1t NaOH = 25000mol NaOH = 2.5e7 meq Alk
duration = 1hour
depth = -1
radius = 50
release_rate = total_release/duration/(π*radius^2*abs(depth))

oae = Forcing(alkalinity_release; 
              parameters = (; start_time = 20minutes, 
                              duration, 
                              release_rate, 
                              radius, 
                              depth)) 

# And then put the model together with `Alk_release` and `Alk_no_precipitation` both receiving
# the release while `Alk_control` evolves normally

model = NonhydrostaticModel(grid;
                            coriolis = FPlane(latitude = 45),
                            boundary_conditions,
                            biogeochemistry,
                            advection = WENO(),
                            tracers = (:T, :S),
                            forcing = (; Alk_release = oae, Alk_no_precipitation = oae))

set!(model; P = 1, NO₃ = 10,
            T = 15, S = 35,
            (replicate.DIC => 2000 for replicate in replicates)...,
            (replicate.Alk => 2300 for replicate in replicates)...,
            (replicate.CaCO₃ => 0.1 for replicate in replicates)...,
            u = (x, y, z) -> randn()/2)

# Then we can setup a simulation with a variable timestep and progress message
# and output the tracers, running for longer than the release itself so the air-sea flux
# has time to respond

simulation = Simulation(model, Δt = 5, stop_time = 4hours)

conjure_time_step_wizard!(simulation; cfl = 0.5)

prog(sim) = 
    @info prettytime(sim)*" in "*prettytime(sim.run_wall_time)*" with Δt = "*prettytime(sim.Δt)

add_callback!(simulation, prog, IterationInterval(100))

# we save the saturation states alongside the tracers so the state of the explicit
# calcium carbonate pool can be inspected

outputs = merge(model.tracers, (; (replicate.Ω => model.auxiliary_fields[replicate.Ω] for replicate in replicates)...))

simulation.output_writers[:tracers] = JLD2Writer(model, outputs;
                                                filename = "oae.jld2",
                                                schedule = AveragedTimeInterval(20minutes),
                                                overwrite_files = true)

# We can also create a `BoundaryConditionOperation` which records the flux through the
# top boundary for each of the DIC fields which we can save (as `qCO₂_control`, ...)

qCO₂ = (; (Symbol(:qCO₂_, replicate.label) => BoundaryConditionOperation(model.tracers[replicate.DIC], :top, model)
           for replicate in replicates)...)

simulation.output_writers[:carbon_flux] = JLD2Writer(model, qCO₂,
                                                     indices = (:, :, grid.Nz),
                                                     filename = "oae_surface_flux.jld2",
                                                     schedule = TimeInterval(20minutes),
                                                     overwrite_files = true)

# and then run the simulation

run!(simulation)

# Finally we load the data and create some plots

using CairoMakie, Statistics

fds = FieldDataset("oae.jld2")
fds_surface = FieldDataset("oae_surface_flux.jld2")

n = Observable(1)

N = length(fds["P"])

# compare the `release` replicate (precipitation active) with the `control`

Δ_Alk = @lift (mean(interior(fds["Alk_release"][$n]), dims=3)[:, :, 1] .-
               mean(interior(fds["Alk_control"][$n]), dims=3)[:, :, 1]).* 15
Δ_qCO₂ = @lift (interior(fds_surface["qCO₂_release"][$n], :, :, 1) .-
                interior(fds_surface["qCO₂_control"][$n], :, :, 1)) .* (12+2*16)*1e-3*day

Δ_CaCO₃ = @lift (mean(interior(fds["CaCO₃_release"][$n]), dims=3)[:, :, 1] .-
                 mean(interior(fds["CaCO₃_control"][$n]), dims=3)[:, :, 1]) .* 15

Δ_Alk_range = maximum(abs, mean(interior(fds["Alk_release"]), dims = 3) .-
                           mean(interior(fds["Alk_control"]), dims = 3)) .* 15
Δ_CaCO₃_range = maximum(abs, mean(interior(fds["CaCO₃_release"]), dims = 3) .-
                             mean(interior(fds["CaCO₃_control"]), dims = 3)) .* 15
Δ_qCO₂_range = maximum(abs, interior(fds_surface["qCO₂_release"]) .-
                            interior(fds_surface["qCO₂_control"])) .* (12+2*16)*1e-3*day

xc = xnodes(fds["P"])
yc = ynodes(fds["P"])

fig = Figure(size=(1500, 500))

ax = Axis(fig[1, 1], aspect = DataAspect())
ax2 = Axis(fig[1, 2], aspect = DataAspect())
ax3 = Axis(fig[1, 3], aspect = DataAspect())

hm = heatmap!(ax, xc, yc, Δ_Alk, colormap = :balance, colorrange = (-1, 1) .* Δ_Alk_range)
hm2 = heatmap!(ax2, xc, yc, Δ_qCO₂, colormap = :balance, colorrange = (-1, 1) .* Δ_qCO₂_range)
hm3 = heatmap!(ax3, xc, yc, Δ_CaCO₃, colormap = :balance, colorrange = (-1, 1) .* Δ_CaCO₃_range)

Colorbar(fig[2, 1], hm, vertical = false, 
                        label = "Alkalinity pertubation (mmol / m²)", 
                        flip_vertical_label = true)
Colorbar(fig[2, 2], hm2, vertical = false, 
                         label = "Surface CO₂ flux difference (gCO₂/m²/day)", 
                         flip_vertical_label = true)
Colorbar(fig[2, 3], hm3, vertical = false,
                         label = "Calcium carbonate difference (mmol / m²)",
                         flip_vertical_label = true)

supertitle = Label(fig[0, :], (@lift prettytime(fds["P"].times[$n])))

CairoMakie.record(fig, "oae.mp4", 1:N, framerate=10) do i
    @info "$i"
    n[] = i
end

# ![](oae.mp4)

fig = Figure(size = (1000, 750));

ax = Axis(fig[1, 1], title = "Surface flux (gC/day)")

surface_area = 500*500
mmol_to_g = 12/1000

mean_flux(name) = mean(interior(fds_surface[name], :, :, 1, :), dims = (1, 2))[1, 1, :]

control_flux = mean_flux("qCO₂_control")
no_precipitation_flux = mean_flux("qCO₂_no_precipitation")
release_flux = mean_flux("qCO₂_release")

times = fds_surface["qCO₂_control"].times ./ hours

lines!(ax, times, control_flux .* mmol_to_g * surface_area * day, label = "Control")
lines!(ax, times, no_precipitation_flux .* mmol_to_g * surface_area * day, label = "Release, no precipitation")
lines!(ax, times, release_flux .* mmol_to_g * surface_area * day, label = "Release, precipitation",
       linestyle = :dash)

Legend(fig[1, 2], ax)

# Cumulative difference in the air-sea CO₂ flux between replicate pairs

ax2 = Axis(fig[2, 1], title = "Cumulative CO₂ difference (gC)")

Δt_save = diff(fds_surface["qCO₂_control"].times)[1]

# a negative flux is carbon entering the ocean, so `a - b` is the extra carbon taken up by `b`
cumulative_uptake(a, b) = cumsum(Δt_save .* (mean_flux(a) .- mean_flux(b)) .* mmol_to_g * surface_area)

lines!(ax2, times, cumulative_uptake("qCO₂_control", "qCO₂_no_precipitation"),
       label = "Extra CO₂ taken up, no precipitation")
lines!(ax2, times, cumulative_uptake("qCO₂_control", "qCO₂_release"),
       label = "Extra CO₂ taken up, precipitation",
       linestyle = :dash)

Legend(fig[2, 2], ax2)

domain_volume = surface_area * 15

# Finally, the fate of the added alkalinity itself: how much has been realised as CO₂ uptake versus
# destroyed by precipitation. The uptake is converted to an alkalinity equivalent using the
# equilibrium efficiency ``\eta = \Delta DIC / \Delta TA``, found by asking what `DIC` equilibrates
# with the atmosphere before and after an addition:

function equilibrium_DIC(Alk; target = 413.0, T = 15.0, S = 35.0)
    lo, hi = 1000.0, 3500.0
    for _ in 1:60
        mid = (lo + hi) / 2
        carbon_chemistry(; DIC = mid, Alk, T, S) < target ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

η = (equilibrium_DIC(2400.0) - equilibrium_DIC(2300.0)) / 100

ax3 = Axis(fig[3, 1], xlabel = "Time (hours)", ylabel = "Alkalinity (meq)",
                      title = "Fate of the added alkalinity (release vs control)")

TA_realised = cumsum(Δt_save .* (mean_flux("qCO₂_control") .- mean_flux("qCO₂_release")) .* surface_area) ./ η
TA_precipitated = [2 * mean(interior(fds["CaCO₃_release"][i]) .- interior(fds["CaCO₃_control"][i])) * domain_volume
                   for i in 1:N]

lines!(ax3, times, TA_realised, label = "Converted to CO₂ uptake")
lines!(ax3, fds["CaCO₃_control"].times./hours, TA_precipitated, label = "Lost to precipitation")

Legend(fig[3, 2], ax3)

fig

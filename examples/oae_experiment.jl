# # Idealised ocean alkalinity enhancement (OAE)
#
# In this example we setup a 3D model driven by a constant wind stress, with 2 instances
# of the carbon chemistry system (DIC and alkalinity). We then add alkalinity to one of
# the instances replicating an OAE release, and compare to the unperturbed counterfactual

# ## Install dependencies
# First we ensure we have the required dependencies installed
# ```julia
# using Pkg
# pkg "add OceanBioME, Oceananigans, CairoMakie"
# ```

# ## Model setup

using Oceananigans, OceanBioME, Oceananigans.Units

using OceanBioME.Models.GasExchangeModel: CarbonDioxideConcentration

grid = RectilinearGrid(GPU(); 
                       size = (64, 64, 8), 
                       extent = (500, 500, 15))

biogeochemistry = LOBSTER(grid;
                          carbonate_system = CarbonateSystem(2))

carbon_chemistry = CarbonChemistry()
CO₂_flux1 = CarbonDioxideGasExchangeBoundaryCondition(; water_concentration = 
                                                            CarbonDioxideConcentration(; carbon_chemistry, 
                                                                                         DIC = :DIC1,
                                                                                         Alk = :Alk1))
CO₂_flux2 = CarbonDioxideGasExchangeBoundaryCondition(; water_concentration = 
                                                            CarbonDioxideConcentration(; carbon_chemistry, 
                                                                                         DIC = :DIC2,
                                                                                         Alk = :Alk2))
wind = FluxBoundaryCondition(-1.2/1026*2e-3*10^2)

boundary_conditions = (; DIC1 = FieldBoundaryConditions(top = CO₂_flux1),
                         DIC2 = FieldBoundaryConditions(top = CO₂_flux2),
                         u    = FieldBoundaryConditions(top = wind))

@inline oae_release(x, y, z, t, params) = 
    ifelse((params.start_time <= t < params.start_time + params.duration) & (z > params.depth) & ((x-250)^2 + (y-250)^2 < params.radius^2), params.release_rate, 0)

total_release = 2.5e7 # 1t NaOH = 25000mol NaOH = 2.5e7 meq Alk
duration = 1hour
depth = -1
radius = 50
release_rate = total_release/duration/(π*radius^2*abs(depth))

oae = Forcing(oae_release; parameters = (; start_time = 1hour, duration, release_rate, radius, depth)) 

model = NonhydrostaticModel(grid;
                            coriolis = FPlane(latitude = 45),
                            boundary_conditions,
                            biogeochemistry,
                            advection = WENO(),
                            tracers = (:T, :S),
                            forcing = (; Alk2 = oae))

set!(model, P = 1, NO₃ = 10, T = 15, S = 35, DIC1 = 2000, Alk1 = 2300, DIC2 = 2000, Alk2 = 2300, u = (x, y, z) -> randn()/2)

simulation = Simulation(model, Δt = 5, stop_time = 5hours)

conjure_time_step_wizard!(simulation)

prog(sim) = @info prettytime(sim) * " in " * prettytime(sim.run_wall_time) * " with Δt = " * prettytime(sim.Δt)

add_callback!(simulation, prog, IterationInterval(100))

simulation.output_writers[:tracers] = JLD2Writer(model, model.tracers;
                                                filename = "oae.jld2",
                                                schedule = AveragedTimeInterval(5minutes),
                                                overwrite_existing = true)

qCO₂1 = BoundaryConditionOperation(model.tracers.DIC1, :top, model)
qCO₂2 = BoundaryConditionOperation(model.tracers.DIC2, :top, model)


simulation.output_writers[:carbon_flux] = JLD2Writer(model, (; qCO₂1, qCO₂2),
                                                     indices = (:, :, grid.Nz),
                                                     filename = "oae_surface_flux.jld2",
                                                     schedule = TimeInterval(5minutes),
                                                     overwrite_existing = true)

run!(simulation)

using CairoMakie, Statistics

fds = FieldDataset("oae.jld2")
fds_surface = FieldDataset("oae_surface_flux.jld2")

n = Observable(1)

N = length(fds["P"])

Δ_Alk = @lift (mean(interior(fds["Alk2"][$n]), dims=3)[:, :, 1] .- mean(interior(fds["Alk1"][$n]), dims=3)[1, 1, :]).* 15
Δ_qCO₂ = @lift (interior(fds_surface["qCO₂2"][$n], :, :, 1) .- interior(fds_surface["qCO₂1"][$n], :, :, 1)) .* (12+2*16)*1e-3*day

Δ_Alk_range = maximum(abs, mean(interior(fds["Alk2"]), dims = 3) .- mean(interior(fds["Alk1"]), dims = 3)) .* 15
Δ_qCO₂_range = maximum(abs, interior(fds_surface["qCO₂2"]) .- interior(fds_surface["qCO₂1"])) .* (12+2*16)*1e-3*day

xc = xnodes(fds["P"])
yc = ynodes(fds["P"])

fig = Figure(size=(1000, 500))

ax = Axis(fig[1, 1], aspect = DataAspect())
ax2 = Axis(fig[1, 2], aspect = DataAspect())

hm = heatmap!(ax, xc, yc, Δ_Alk, colormap = :balance, colorrange = (-1, 1) .* Δ_Alk_range)
hm2 = heatmap!(ax2, xc, yc, Δ_qCO₂, colormap = :balance, colorrange = (-1, 1) .* Δ_qCO₂_range)

Colorbar(fig[2, 1], hm, vertical = false, label = "Alkalinity pertubation (mmol / m³)", flip_vertical_label = true)
Colorbar(fig[2, 2], hm2, vertical = false, label = "Surface CO₂ flux difference (gCO₂/m²/day)", flip_vertical_label = true)

supertitle = Label(fig[0, :], (@lift prettytime(fds["P"].times[$n])))

CairoMakie.record(fig, "oae.mp4", 1:N, framerate=10) do i
    @info "$i"
    n[] = i
end

# ![](oae.mp4)

fig = Figure();

ax = Axis(fig[1, 1], title = "Surface flux (gC/day)")

surface_area = 500*500
mmol_to_g = 12/1000

lines!(ax, fds_surface["qCO₂1"].times./hours, mean(interior(fds_surface["qCO₂1"], :, :, 1, :), dims = (1, 2))[1, 1, :].* mmol_to_g * surface_area * day)
lines!(ax, fds_surface["qCO₂2"].times./hours, mean(interior(fds_surface["qCO₂2"], :, :, 1, :), dims = (1, 2))[1, 1, :].* mmol_to_g * surface_area * day)

ax2 = Axis(fig[2, 1], xlabel = "Time (hours)", title = "Cumulative difference (gC)")

Δt = diff(fds_surface["qCO₂1"].times./hours)[1]

change =  mean(interior(fds_surface["qCO₂1"], :, :, 1, :), dims = (1, 2))[1, 1, :] .- mean(interior(fds_surface["qCO₂2"], :, :, 1, :), dims = (1, 2))[1, 1, :]
change .*= mmol_to_g * surface_area

Δt_save = diff(fds_surface["qCO₂1"].times)[1]

cumulative_change = cumsum(Δt_save .* change)

lines!(ax2, fds_surface["qCO₂1"].times./hours, cumulative_change)

fig
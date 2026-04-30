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

total_release = 2.5e8 # 1t NaOH
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

set!(model, P = 1, NO₃ = 10, T = 15, S = 35, DIC1 = 1035.5, Alk1 = 1100.065, DIC2 = 1035.5, Alk2 = 1100.065, u = (x, y, z) -> randn()/2)

simulation = Simulation(model, Δt = 5, stop_time = 1day)

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

using CairoMakie

fds = FieldDataset("oae.jld2")
fds_surface = FieldDataset("oae_surface_flux.jld2")

n = Observable(1)

N = length(fds["P"])

P_plt = @lift view(fds["P"][$n], :, :, grid.Nz)
Δ_Alk = @lift interior(fds["Alk2"][$n], :, :, grid.Nz) .- interior(fds["Alk1"][$n], :, :, grid.Nz)
Δ_qCO₂ = @lift (interior(fds_surface["qCO₂2"][$n], :, :, 1) .- interior(fds_surface["qCO₂1"][$n], :, :, 1)) .* (12+2*16)*1e-3*day
Δ_DIC = @lift interior(fds["DIC2"][$n], :, :, Alk1.grid.Nz) .- interior(fds["DIC1"][$n], :, :, Alk1.grid.Nz)

P_range = (minimum(fds["P"]), maximum(fds["P"]))
Δ_Alk_range = maximum(map(n->maximum(abs, fds["Alk2"][n] - fds["Alk1"][n]), 1:N))
Δ_qCO₂_range = maximum(map(n->maximum(abs, fds_surface["qCO₂2"][n] -fds_surface["qCO₂1"][n]), 1:N) .* (12+2*16)*1e-3*day)
Δ_DIC_range = maximum(map(n->maximum(abs, fds["DIC2"][n] - fds["DIC1"][n]), 1:N))

xc = xnodes(fds["P"])
yc = ynodes(fds["P"])

fig = Figure(size=(1000, 1150))

ax = Axis(fig[2, 1], aspect = DataAspect())
ax2 = Axis(fig[2, 2], aspect = DataAspect())
ax3 = Axis(fig[4, 1], aspect = DataAspect())
ax4 = Axis(fig[4, 2], aspect = DataAspect())

hm1 = heatmap!(ax, P_plt, colorrange = P_range)
hm2 = heatmap!(ax2, xc, yc, Δ_Alk, colormap = :balance, colorrange = (-1, 1) .* Δ_Alk_range)
hm3 = heatmap!(ax3, xc, yc, Δ_qCO₂, colormap = :balance, colorrange = (-1, 1) .* Δ_qCO₂_range)
hm4 = heatmap!(ax4, xc, yc, Δ_DIC, colormap = :balance, colorrange = (-1, 1) .* Δ_DIC_range)

Colorbar(fig[1, 1], hm1, vertical = false, label = "Phytoplankton (mmol N/m³)")
Colorbar(fig[1, 2], hm2, vertical = false, label = "Alkalinity pertubation (mmol / m³)")
Colorbar(fig[3, 1], hm3, vertical = false, label = "Surface CO₂ flux difference (gCO₂/m²/day)")
Colorbar(fig[3, 2], hm4, vertical = false, label = "DIC pertubation (mmol C / m³)")

supertitle = Label(fig[0, :], (@lift prettytime(fds["P"].times[$n])))

CairoMakie.record(fig, "oae.mp4", 1:N, framerate=10) do i
    @info "$i"
    n[] = i
end

fig = Figure()

ax = Axis(fig[1, 1], ylabel = "Surface CO₂ flux (gCO₂/m²/day)")

lines!(ax, fds_surface["qCO₂1"].times./day, map(n->mean(interior(fds_surface["qCO₂1"][n], :, :, 1)).* (12+2*16)*1e-3*day, 1:N))
lines!(ax, fds_surface["qCO₂2"].times./day, map(n->mean(interior(fds_surface["qCO₂2"][n], :, :, 1)).* (12+2*16)*1e-3*day, 1:N))

ax2 = Axis(fig[2, 1], xlabel = "Time (hours)", ylabel = "Cumulative CO₂ flux difference (gCO₂/m²)")

Δt = diff(fds_surface["qCO₂1"].times./day)[1]

lines!(ax2, fds_surface["qCO₂1"].times./day, Δt .* cumsum(map(n->(mean(interior(fds_surface["qCO₂2"][n], :, :, 1)) - mean(interior(fds_surface["qCO₂1"][n], :, :, 1))).* (12+2*16)*1e-3*day, 1:N)))
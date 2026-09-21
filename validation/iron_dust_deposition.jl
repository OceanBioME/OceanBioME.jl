# Column validation of IronDustDeposition + SimpleIron
#
# Runs two 1-D column simulations using the MITgcmDIC preset (SimpleIron + ImplicitProductivity):
#   1. No dust — iron is only recycled and scavenged
#   2. With dust — constant aeolian dust deposition supplies iron to the surface ocean
#
# The comparison shows how dust-derived iron relieves iron limitation and sustains productivity.

using OceanBioME, Oceananigans, Printf
using Oceananigans.Fields: FunctionField, ConstantField
using Oceananigans.Units

const year = years = 365days

# ── physical forcing ──────────────────────────────────────────────────────────

@inline PAR⁰(t) = 60 * (1 - cos((t + 15days) * 2π / year)) * (1 / (1 + 0.2 * exp(-((mod(t, year) - 200days) / 50days)^2))) + 2

@inline H(t, t₀, t₁) = ifelse(t₀ < t < t₁, 1.0, 0.0)
@inline fmld1(t) = H(t, 50days, year) * (1 / (1 + exp(-(t - 100days) / 5days))) * (1 / (1 + exp((t - 330days) / 25days)))
@inline MLD(t) = -(10 + 340 * (1 - fmld1(year - eps(year)) * exp(-mod(t, year) / 25days) - fmld1(mod(t, year))))
@inline κₜ(z, t) = 1e-2 * (1 + tanh((z - MLD(t)) / 10)) / 2 + 1e-4

# ── grid ──────────────────────────────────────────────────────────────────────

grid = RectilinearGrid(topology = (Flat, Flat, Bounded), size = (50, ), extent = (500, ))

# ── run a single column ──────────────────────────────────────────────────────

function run_column(grid; dust_deposition = nothing, stop_time = 3years, label = "column")
    clock = Clock(; time = 0.0)
    κ_field = FunctionField{Center, Center, Center}(κₜ, grid; clock)

    forcing = if isnothing(dust_deposition)
        NamedTuple()
    else
        (; Fe = IronDustDepositionForcing(dust_deposition))
    end

    biogeochemistry = MITgcmDIC(grid;
                                surface_PAR = PAR⁰,
                                oxygen = nothing)

    CO₂_flux = CarbonDioxideGasExchangeBoundaryCondition()

    T = ConstantField(10.0)
    S = ConstantField(35.0)

    model = HydrostaticFreeSurfaceModel(grid;
                velocities = PrescribedVelocityFields(),
                tracer_advection = nothing,
                momentum_advection = nothing,
                buoyancy = nothing,
                clock,
                closure = ScalarDiffusivity(VerticallyImplicitTimeDiscretization(), κ = κ_field),
                biogeochemistry,
                boundary_conditions = (DIC = FieldBoundaryConditions(top = CO₂_flux),),
                auxiliary_fields = (; T, S),
                forcing)

    set!(model, PO₄ = 1.5, Fe = 5e-5, DOP = 0.0, POP = 0.0, DIC = 2200.0, Alk = 2400.0)

    simulation = Simulation(model, Δt = 10minutes, stop_time = stop_time)

    progress_message(sim) = @printf("  [%s] Iteration: %04d, time: %s, Δt: %s, wall time: %s\n",
                                    label, iteration(sim), prettytime(sim),
                                    prettytime(sim.Δt), prettytime(sim.run_wall_time))

    add_callback!(simulation, progress_message, TimeInterval(30day))

    filename = label
    simulation.output_writers[:tracers] = JLD2Writer(model, model.tracers,
                                                     filename = "$filename.jld2",
                                                     schedule = TimeInterval(1day),
                                                     overwrite_files = true)

    @info "Running $label..."
    run!(simulation)

    return filename
end

# ── run both cases ────────────────────────────────────────────────────────────

# typical global-mean dust deposition: ~2 g/m²/year ≈ 6.3e-11 kg/m²/s
dust_flux = 2e-3 / year  # kg/m²/s

file_nodust = run_column(grid; dust_deposition = nothing,   label = "iron_nodust")
file_dust   = run_column(grid; dust_deposition = dust_flux, label = "iron_dust")

# ── load and plot ─────────────────────────────────────────────────────────────

using CairoMakie

function load_timeseries(filename, tracer)
    return FieldTimeSeries(filename * ".jld2", tracer)
end

Fe_nodust  = load_timeseries(file_nodust, "Fe")
Fe_dust    = load_timeseries(file_dust,   "Fe")
PO₄_nodust = load_timeseries(file_nodust, "PO₄")
PO₄_dust   = load_timeseries(file_dust,   "PO₄")
DIC_nodust = load_timeseries(file_nodust, "DIC")
DIC_dust   = load_timeseries(file_dust,   "DIC")

x, y, z = nodes(Fe_nodust)
times = Fe_nodust.times

fig = Figure(size = (1400, 1200), fontsize = 18)

start_day = 1
end_day   = length(times)
axis_kwargs = (xlabel = "Time (days)", ylabel = "z (m)",
               limits = ((times[start_day], times[end_day]) ./ days, (-300, 0)))

# ── shared color ranges for fair comparison ──────────────────────────────────
Fe_all = vcat(vec(interior(Fe_nodust)), vec(interior(Fe_dust)))
PO₄_all = vcat(vec(interior(PO₄_nodust)), vec(interior(PO₄_dust)))
DIC_all = vcat(vec(interior(DIC_nodust)), vec(interior(DIC_dust)))

Fe_range  = (minimum(Fe_all),  maximum(Fe_all))
PO₄_range = (minimum(PO₄_all), maximum(PO₄_all))
DIC_range = (minimum(DIC_all), maximum(DIC_all))

# ── Fe heatmaps ──────────────────────────────────────────────────────────────
ax1 = Axis(fig[1, 1]; title = "Fe — no dust (mmol/m³)", axis_kwargs...)
hm1 = heatmap!(ax1, times[start_day:end_day] ./ days, z,
               interior(Fe_nodust, 1, 1, :, start_day:end_day)', colormap = :viridis, colorrange = Fe_range)
lines!(ax1, times[start_day:end_day] ./ days, t -> MLD(t * day), color = :white, linestyle = :dash)
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[1, 3]; title = "Fe — with dust (mmol/m³)", axis_kwargs...)
hm2 = heatmap!(ax2, times[start_day:end_day] ./ days, z,
               interior(Fe_dust, 1, 1, :, start_day:end_day)', colormap = :viridis, colorrange = Fe_range)
lines!(ax2, times[start_day:end_day] ./ days, t -> MLD(t * day), color = :white, linestyle = :dash)
Colorbar(fig[1, 4], hm2)

# ── PO₄ heatmaps ─────────────────────────────────────────────────────────────
ax3 = Axis(fig[2, 1]; title = "PO₄ — no dust (mmol/m³)", axis_kwargs...)
hm3 = heatmap!(ax3, times[start_day:end_day] ./ days, z,
               interior(PO₄_nodust, 1, 1, :, start_day:end_day)', colormap = :batlow, colorrange = PO₄_range)
lines!(ax3, times[start_day:end_day] ./ days, t -> MLD(t * day), color = :white, linestyle = :dash)
Colorbar(fig[2, 2], hm3)

ax4 = Axis(fig[2, 3]; title = "PO₄ — with dust (mmol/m³)", axis_kwargs...)
hm4 = heatmap!(ax4, times[start_day:end_day] ./ days, z,
               interior(PO₄_dust, 1, 1, :, start_day:end_day)', colormap = :batlow, colorrange = PO₄_range)
lines!(ax4, times[start_day:end_day] ./ days, t -> MLD(t * day), color = :white, linestyle = :dash)
Colorbar(fig[2, 4], hm4)

# ── DIC heatmaps ──────────────────────────────────────────────────────────────
ax5 = Axis(fig[3, 1]; title = "DIC — no dust (mmol/m³)", axis_kwargs...)
hm5 = heatmap!(ax5, times[start_day:end_day] ./ days, z,
               interior(DIC_nodust, 1, 1, :, start_day:end_day)', colormap = :tempo, colorrange = DIC_range)
lines!(ax5, times[start_day:end_day] ./ days, t -> MLD(t * day), color = :white, linestyle = :dash)
Colorbar(fig[3, 2], hm5)

ax6 = Axis(fig[3, 3]; title = "DIC — with dust (mmol/m³)", axis_kwargs...)
hm6 = heatmap!(ax6, times[start_day:end_day] ./ days, z,
               interior(DIC_dust, 1, 1, :, start_day:end_day)', colormap = :tempo, colorrange = DIC_range)
lines!(ax6, times[start_day:end_day] ./ days, t -> MLD(t * day), color = :white, linestyle = :dash)
Colorbar(fig[3, 4], hm6)

# ── surface time series comparison ────────────────────────────────────────────
ax7 = Axis(fig[4, 1:2]; xlabel = "Time (days)", ylabel = "Surface Fe (mmol/m³)",
           title = "Surface dissolved iron")
lines!(ax7, times ./ days, interior(Fe_nodust, 1, 1, grid.Nz, :), label = "No dust", linewidth = 2)
lines!(ax7, times ./ days, interior(Fe_dust,   1, 1, grid.Nz, :), label = "With dust", linewidth = 2)
axislegend(ax7, position = :rt)

ax8 = Axis(fig[4, 3:4]; xlabel = "Time (days)", ylabel = "Surface PO₄ (mmol/m³)",
           title = "Surface phosphate")
lines!(ax8, times ./ days, interior(PO₄_nodust, 1, 1, grid.Nz, :), label = "No dust", linewidth = 2)
lines!(ax8, times ./ days, interior(PO₄_dust,   1, 1, grid.Nz, :), label = "With dust", linewidth = 2)
axislegend(ax8, position = :rt)

save("validation/iron_dust_deposition_validation.png", fig, px_per_unit = 2)

fig

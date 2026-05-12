# # CaCO₃ Precipitation, Sinking, and Dissolution in a 1-D column model

# This example demonstrates the `SimpleCaCO3Precipitation` biogeochemical model in a single column model. Temperature, salinity, DIC, and alkalinity are initialized with idealized profiles with a 100m mixed layer and exponentially decaying profiles below. The alkalinity is then perturbed by 200 mmol m⁻³ of alkalinity throughout the 100 m mixed layer. Note that since this is a column model, the added alkalinity does not decrease due to horizontal mixing. All particles that form through secondary precipitation are assumed to be aragonite, which is more soluble than calcite. Note, however, that a large dissolution rate is used for illustration/testing purposes.

# ## Install dependencies
# ```julia
# using Pkg
# pkg"add OceanBioME, Oceananigans, CairoMakie"
# ```

using Printf: @printf, @sprintf
using Statistics: quantile
using OceanBioME, Oceananigans, Oceananigans.Units
using Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields
using Oceananigans.Advection: UpwindBiased
using Oceananigans.TurbulenceClosures: VerticallyImplicitTimeDiscretization
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry, KSP_aragonite

# Use grid stretching with higher resolution near the surface
H = 4000.0   # water depth (m)
Nz = 200         # gridpoints
z_stretch = 1.28 # stretching factor
z_faces = Float64[-H * (1 - (k / Nz)^z_stretch) for k in 0:Nz]

# Create a grid
grid = RectilinearGrid(CPU();
                       size     = Nz,
                       z        = z_faces,
                       topology = (Flat, Flat, Bounded))

# Depth-dependent vertical diffusivity
d_ML  = 100.0    # mixed-layer depth (m)
κ_ml   = 1e-2   # [m² s⁻¹]  mixed-layer turbulent diffusivity
κ_deep = 1e-5   # [m² s⁻¹]  interior diapycnal diffusivity
@inline κ_func(z, _) = κ_ml * (1 + tanh((z + d_ML) / 10)) / 2 + κ_deep
closure = VerticalScalarDiffusivity(VerticallyImplicitTimeDiscretization(), ν = κ_func, κ = κ_func)

k_d_column     = 0.1 / day
w_CaCO₃_column = 20 / day

bgc = SimpleCaCO3Precipitation(;
    grid,
    carbon_chemistry = CarbonChemistry(; calcite_solubility = KSP_aragonite()),
    sinking_speed             = w_CaCO₃_column,
    dissolution_rate_constant = k_d_column)

# Set the air-sea CO₂ flux boundary condition
CO₂_bc  = CarbonDioxideGasExchangeBoundaryCondition(; wind_speed = 10.0)  # u₁₀ = 10 m s⁻¹
DIC_bcs = FieldBoundaryConditions(top = CO₂_bc)

# Create the model
model = HydrostaticFreeSurfaceModel(grid;
    biogeochemistry    = bgc,
    buoyancy           = nothing,  
    coriolis           = nothing,
    closure,
    tracer_advection   = UpwindBiased(),
    momentum_advection = nothing,
    boundary_conditions = (DIC = DIC_bcs,))


# Set the initial conditions                       
decay_scale = 500.0 # exponential decay scale in the ocean interior
# Temperature [°C]
T_ML, T_deep = 22.0, 2.0
# Salinity [psu]
S_ML, S_deep = 36.4, 34.85
# DIC [mmol m⁻³]
DIC_ML, DIC_deep = 2020.0, 2330.0
# Background total alkalinity background [meq m⁻³] (extra alkalinity added to mixed layer below)
Alk_ML, Alk_deep = 2305.0, 2433.0

@inline T_initial(z)   = (d = clamp(-z, 0.0, H); d <= d_ML ? T_ML : T_deep + (T_ML - T_deep) * exp(-(d - d_ML) / decay_scale))
@inline S_initial(z)   = (d = clamp(-z, 0.0, H); d <= d_ML ? S_ML : S_deep + (S_ML - S_deep) * exp(-(d - d_ML) / decay_scale))
@inline DIC_initial(z) = (d = clamp(-z, 0.0, H); d <= d_ML ? DIC_ML : DIC_deep + (DIC_ML - DIC_deep) * exp(-(d - d_ML) / decay_scale))

# Alkalinity addition: +200 mmol m⁻³ of alkalinity uniformly in the mixed layer.
# This corresponds to 50,000 mol added to a 50 m × 50 m patch of ocean.
ΔAlk_ml = 200.0   # [mmol m⁻³]
@inline Alk_initial(z) = (d = clamp(-z, 0.0, H); (d <= d_ML ? Alk_ML : Alk_deep + (Alk_ML - Alk_deep) * exp(-(d - d_ML) / decay_scale)) + ifelse(z > -d_ML, ΔAlk_ml, 0.0))

set!(model;
     T   = T_initial,
     S   = S_initial,
     DIC = DIC_initial,
     Alk = Alk_initial)

simulation = Simulation(model; Δt = 5minutes, stop_time = 360days)
wizard = TimeStepWizard(cfl = 0.5, max_Δt = 30minutes)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(20))

progress_message(sim) = @printf(
    "day %6.1f / 360,  Δt = %s,  wall time = %s\n",
    sim.model.clock.time / day,
    prettytime(sim.Δt),
    prettytime(sim.run_wall_time))

simulation.callbacks[:progress] = Callback(progress_message, TimeInterval(10days))

# Output: tracers plus saturation state Ω (aragonite solubility product in `CarbonChemistry` above)
output_filename = "caco3_precipitation_column.jld2"

outputs = merge(model.tracers, biogeochemical_auxiliary_fields(model.biogeochemistry))
simulation.output_writers[:column] = JLD2Writer(
    model, outputs;
    filename           = output_filename,
    schedule           = TimeInterval(12hours),
    overwrite_existing = true)

run!(simulation)

using CairoMakie

DIC_ts   = FieldTimeSeries(output_filename, "DIC")
Alk_ts   = FieldTimeSeries(output_filename, "Alk")
Ω_ts     = FieldTimeSeries(output_filename, "Ω")
CaCO₃_ts = FieldTimeSeries(output_filename, "CaCO₃")

times = DIC_ts.times
z     = znodes(grid, Center())
t_d   = times ./ day   # convert to days

# For a 1×1×Nz column the interior is (1, 1, Nz, Nt).
# Transposing the z–t slice gives (Nt, Nz) for heatmap!(ax, times, depths, data).
DIC_data   = interior(DIC_ts,   1, 1, :, :)'
Alk_data   = interior(Alk_ts,   1, 1, :, :)'
Ω_data     = interior(Ω_ts,     1, 1, :, :)'
CaCO₃_data = interior(CaCO₃_ts, 1, 1, :, :)'

# Dissolution rate [mmol C m⁻³ day⁻¹], same law as `SimpleCaCO3Precipitation` with m = 1
# (R_diss = k_d * max(0, 1-Ω) * CaCO₃ in model units mmol C m⁻³ s⁻¹).
R_diss_day = k_d_column .* max.(0, 1 .- Ω_data) .* CaCO₃_data .* day

fig = Figure(size = (1200, 1100), fontsize = 14)

ax_kwargs = (xlabel = "Time (days)", ylabel = "Depth (m)",
             limits = ((0, maximum(t_d)), (-H, 0)))

# Top row: DIC, Alk — bottom row: Ω (aragonite), particulate CaCO₃
ax_DIC   = Axis(fig[1, 1]; title = "DIC (mmol C m⁻³)",              ax_kwargs...)
ax_Alk   = Axis(fig[1, 3]; title = "Total alkalinity (meq m⁻³)",   ax_kwargs...)
ax_Ω     = Axis(fig[2, 1]; title = "Saturation state Ω (aragonite)", ax_kwargs...)
ax_CaCO₃ = Axis(fig[2, 3]; title = "Particulate CaCO₃ (mmol m⁻³)",  ax_kwargs...)

hm_DIC   = heatmap!(ax_DIC,   t_d, z, DIC_data;   colormap = :viridis)
hm_Alk   = heatmap!(ax_Alk,   t_d, z, Alk_data;   colormap = :balance)
hm_Ω     = heatmap!(ax_Ω,     t_d, z, Ω_data;     colormap = :RdBu, colorrange = (0.6, 1.4))
hm_CaCO₃ = heatmap!(ax_CaCO₃, t_d, z, CaCO₃_data;
                    colormap = :YlOrRd_9,
                    colorrange = (0, max(1e-6, quantile(vec(CaCO₃_data), 0.995))))

ax_R = Axis(fig[3, 1:3];
            title = @sprintf("CaCO₃ dissolution rate (mmol C m⁻³ day⁻¹); k_d = %.3g/day, w_sink = %.3g m/day",
                             k_d_column * day, w_CaCO₃_column * day),
            ax_kwargs...)
R_hi = max(1e-12, quantile(vec(R_diss_day), 0.995))
hm_R = heatmap!(ax_R, t_d, z, R_diss_day; colormap = :magma, colorrange = (0, R_hi))

Colorbar(fig[1, 2], hm_DIC;   label = "mmol m⁻³", width = 15)
Colorbar(fig[1, 4], hm_Alk;   label = "meq m⁻³",  width = 15)
Colorbar(fig[2, 2], hm_Ω;     label = "Ω",        width = 15)
Colorbar(fig[2, 4], hm_CaCO₃; label = "mmol m⁻³", width = 15)
Colorbar(fig[3, 4], hm_R;     label = "mmol m⁻³ day⁻¹", width = 15)

figure_filename = "caco3_precipitation_example.png"
save(figure_filename, fig)
@info "Figure saved → $figure_filename"

rm(output_filename; force = true)
@info "JLD2 output deleted"

fig

# Validation of SimpleCaCO3Precipitation against shipboard OAE experiments
#
# Reference:
#   Hashim, M. S. et al. (2025): Mineral formation during shipboard ocean
#   alkalinity enhancement experiments in the North Atlantic.
#   Biogeosciences, 22, 7149–7165. https://doi.org/10.5194/bg-22-7149-2025
#
# Experiment design (Table 1):
#   Sargasso Sea seawater: TA₀ = 2547, DIC₀ = 2082 µmol kg⁻¹, T = 27 °C, S = 36 psu
#   NaOH added to sealed 3 L bags (closed to atmosphere, no DIC change at t = 0).
#   A (control): +0 µmol kg⁻¹
#   B:           +500 µmol kg⁻¹
#   C:           +1000 µmol kg⁻¹
#   D:           +2000 µmol kg⁻¹
#   Run for ~5 days.
#
# Empirical precipitation rate law (Eq. 5 / Fig. 9):
#   log₁₀(r) = n × log₁₀(ΩA − 1) + log₁₀(k)
#   n = 2.15 ± 0.40
#   k = 0.20 ± 0.10 µmol h⁻¹  (for the 3 L bag, unseeded)
#
# Unit conversion for SimpleCaCO3Precipitation (mmol C m⁻³ s⁻¹):
#   Since 1 µmol L⁻¹ = 1 mmol m⁻³, the volumetric rate constant is simply:
#   k_vol = k [µmol h⁻¹] / V_bag [L] / 3600 [s h⁻¹]
#         = 0.20 / 3.0 / 3600 ≈ 1.85 × 10⁻⁵  mmol m⁻³ s⁻¹
#
# Known limitations of the model relative to the paper:
#   1. The model has no induction period. The paper finds induction periods of
#      0.5 – 8 h depending on initial Ω.
#   2. The model has no kinetic inhibition (Mg²⁺, phosphate, organics). The paper
#      observes Ω levelling off at ~3 rather than reaching equilibrium.

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using OceanBioME, Oceananigans, Oceananigans.Units
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry, KSP_aragonite
using JLD2
using Printf

# ---------------------------------------------------------------------------
# Rate parameters  (Hashim et al. 2025, converted to model units)
# ---------------------------------------------------------------------------

const k_paper = 0.20    # µmol h⁻¹  (extensive rate for one 3 L bag)
const V_bag   = 3.0     # L
const n_paper = 2.15    # reaction order

# 1 µmol L⁻¹ = 1 mmol m⁻³, so the volumetric rate constant simplifies to:
const k_vol = k_paper / V_bag / 3600   # mmol m⁻³ s⁻¹

# ---------------------------------------------------------------------------
# Experimental conditions  (Table 1)
# µmol kg⁻¹ used directly as mmol m⁻³ (ρ ≈ 1025 kg m⁻³ → ~2.5 % error)
# ---------------------------------------------------------------------------

const T_exp  = 27.0    # °C
const S_exp  = 36.0    # psu
const TA_bg  = 2547.0  # mmol m⁻³
const DIC_bg = 2082.0  # mmol m⁻³

# Initial particulate CaCO₃ from reported PIC = 0.13 µmol L⁻¹ (Table 1)
const CaCO₃_bg = 0.13  # mmol m⁻³

experiments = (
    (name = "A_control", label = "Control (+0)",   ΔTA =    0.0),
    (name = "B_500",     label = "Exp B (+500)",   ΔTA =  500.0),
    (name = "C_1000",    label = "Exp C (+1000)",  ΔTA = 1000.0),
    (name = "D_2000",    label = "Exp D (+2000)",  ΔTA = 2000.0),
)

# ---------------------------------------------------------------------------
# Run one experiment as an OceanBioME box model
# ---------------------------------------------------------------------------

function run_experiment(exp; stop_time = 5days, Δt = 30seconds)

    grid = BoxModelGrid()

    bgc = SimpleCaCO3Precipitation(;
        grid,
        precipitation_rate_constant = k_vol,
        precipitation_order         = n_paper,
        dissolution_rate_constant   = 0.0,  # sealed bags: precipitate stays in
        sinking_speed               = 0.0,  # 0-D: sinking is meaningless
        carbon_chemistry            = CarbonChemistry(; calcite_solubility = KSP_aragonite()))

    model = BoxModel(; biogeochemistry = bgc, clock = Clock(time = 0.0))

    set!(model; DIC = DIC_bg, Alk = TA_bg + exp.ΔTA, CaCO₃ = CaCO₃_bg, T = T_exp, S = S_exp)

    filename = "hashim2025_$(exp.name).jld2"
    Ω_filename = "hashim2025_$(exp.name)_Omega.jld2"

    simulation = Simulation(model; Δt, stop_time)

    # Tracer timeseries
    simulation.output_writers[:fields] = JLD2Writer(model, model.fields;
                                                     filename,
                                                     schedule = TimeInterval(1hours),
                                                     overwrite_existing = true)

    # Saturation state Ω (updated each timestep by update_biogeochemical_state!)
    Ω_field = model.biogeochemistry.underlying_biogeochemistry.calcite_saturation
    simulation.output_writers[:saturation] = JLD2Writer(model, (; Ω = Ω_field);
                                                         filename = Ω_filename,
                                                         schedule = TimeInterval(1hours),
                                                         overwrite_existing = true)

    progress_interval = round(Int, 12hours / Δt)
    simulation.callbacks[:progress] = Callback(
        sim -> @printf("[%s]  t = %-14s  TA = %7.1f  DIC = %7.1f  Ω = %.3f\n",
                       exp.label,
                       prettytime(time(sim)),
                       sim.model.fields.Alk[1, 1, 1],
                       sim.model.fields.DIC[1, 1, 1],
                       Ω_field[1, 1, 1]),
        IterationInterval(progress_interval))

    @info "Running $(exp.label)…"
    run!(simulation)

    return (; filename, Ω_filename)
end

files = [run_experiment(exp) for exp in experiments]

# ---------------------------------------------------------------------------
# Load timeseries
# ---------------------------------------------------------------------------

function load_ts(f)
    t     = FieldTimeSeries(f.filename, "DIC").times
    DIC   = FieldTimeSeries(f.filename, "DIC")[1, 1, 1, :]
    Alk   = FieldTimeSeries(f.filename, "Alk")[1, 1, 1, :]
    CaCO₃ = FieldTimeSeries(f.filename, "CaCO₃")[1, 1, 1, :]
    Ω     = FieldTimeSeries(f.Ω_filename, "Ω")[1, 1, 1, :]
    return (; t, DIC, Alk, CaCO₃, Ω)
end

data = [load_ts(f) for f in files]

# ---------------------------------------------------------------------------
# Print summary table  (compare to paper's Fig. 6A / 10B)
# ---------------------------------------------------------------------------

@info "\n  Experiment   │  TA initial  │  TA final  │  ΔTA removed  │  Ω initial"
@info "  ─────────────┼──────────────┼────────────┼───────────────┼───────────"
for (exp, ts) in zip(experiments, data)
    @printf("  %-13s│  %9.1f   │  %8.1f  │  %11.1f  │  %.2f\n",
            exp.label, ts.Alk[1], ts.Alk[end], ts.Alk[1] - ts.Alk[end], ts.Ω[1])
end

# ---------------------------------------------------------------------------
# Plots  (mirrors Figs 6A, 6B, 7A, and 8 from the paper)
# ---------------------------------------------------------------------------

using CairoMakie

colors = [:black, :steelblue, :darkorange, :crimson]
h_per_s = 1 / 3600

# Unit conversion: mmol CaCO₃ m⁻³ → mg kg⁻¹
# CaCO₃ [mg/kg] = CaCO₃ [mmol/m³] × M_CaCO₃ [g/mol] / ρ_sw [kg/m³]
const M_CaCO₃ = 100.09   # g/mol
const ρ_sw    = 1025.0   # kg/m³
const mmolm3_to_mgkg = M_CaCO₃ / ρ_sw

fig = Figure(size = (1000, 900), fontsize = 13)

ax_TA    = Axis(fig[1, 1]; xlabel = "Time (h)", ylabel = "TA  (mmol m⁻³)",
                title = "Total Alkalinity  [cf. Fig. 6a]")
ax_DIC   = Axis(fig[1, 2]; xlabel = "Time (h)", ylabel = "DIC  (mmol m⁻³)",
                title = "Dissolved Inorganic Carbon  [cf. Fig. 6b]")
ax_Ω     = Axis(fig[2, 1]; xlabel = "Time (h)", ylabel = "Ω_aragonite  (–)",
                title = "Aragonite Saturation State  [cf. Fig. 7a]")
ax_CaCO₃ = Axis(fig[2, 2]; xlabel = "Time (h)", ylabel = "PIC  (mg kg⁻¹)",
                title = "Particulate Inorganic Carbon  [cf. Fig. 8]")

for (exp, ts, col) in zip(experiments, data, colors)
    t_h = ts.t .* h_per_s

    lines!(ax_TA,    t_h, ts.Alk;   color = col, label = exp.label, linewidth = 2)
    lines!(ax_DIC,   t_h, ts.DIC;   color = col, linewidth = 2)
    lines!(ax_Ω,     t_h, ts.Ω;    color = col, linewidth = 2)
    lines!(ax_CaCO₃, t_h, ts.CaCO₃ .* mmolm3_to_mgkg; color = col, linewidth = 2)
end

# Reference lines
hlines!(ax_TA, [TA_bg]; linestyle = :dash, color = :grey55, label = "Background TA")
hlines!(ax_Ω,  [1.0];   linestyle = :dash, color = :grey55, label = "Ω = 1 (equilibrium)")

for ax in [ax_TA, ax_DIC, ax_Ω, ax_CaCO₃]
    xlims!(ax, 0, 5 * 24)
end

Legend(fig[3, :], ax_TA; orientation = :horizontal, framevisible = false, tellwidth = false)

save("hashim2025_comparison.png", fig)
@info "Figure saved → hashim2025_comparison.png"

# ---------------------------------------------------------------------------
# Clean up temporary JLD2 output files
# ---------------------------------------------------------------------------

for f in files
    rm(f.filename;   force = true)
    rm(f.Ω_filename; force = true)
end
@info "JLD2 output files deleted"

fig

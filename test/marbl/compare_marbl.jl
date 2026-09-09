#####
##### MARBL numerical comparison — the SINGLE consolidated Fortran-comparison script for the whole
##### MARBL-in-NPD model. Run in the isolated env that carries NCDatasets, --check-bounds=yes:
#####   julia --check-bounds=yes --project=<cmpenv> marbl-info3/draft/compare_marbl.jl
#####
##### This replaces the per-phase comparison files (compare_marbl_{Jtend,realmethod,phase5,phase6,phase7,
##### cocco}.jl, archived). It feeds MARBL's exact 60-level column state (the shipped
##### `call_compute_subroutines` baseline) into OUR rate helpers / real tendency methods and matches the
##### standalone MARBL Fortran driver. Sections:
#####   1. per-RATE terms + grazing network (§4/§6/§8/§9/§10/§12/§13/§14) — the phase-3/4 milestone,
#####      incl. the photoC 6-sub-column reconstruction (PCref/alphaPI/εT/thetaN_max/Geider) + neg control
#####   2. REAL METHOD  bgc(Val:X)  vs MARBL J_<tracer>   — autotroph quotas + DOM, machine precision
#####   3. ASSEMBLED  J_<tracer>  reconstructed from MARBL rate diagnostics (validates photoacc synthesis)
#####   4. DOM / POM / PFe rate terms (§15.1/§15.2/§16.1) — production, remin, refractory split, biological PFe
#####   5. carbonate & calcite (§9/§16) — formation/production, J_spCaCO3, carbonate cross-solver
#####   6. O₂ / N redox (§18/§19) — NITRIF, DENITRIF, O2_PROD/CONS, J_O2/J_NO3/J_NH4
#####   8. iron cycle (§17) — speciation / scavenging / ligand budget
#####   9. sediments (§16.5) — pocToSed/ponToSed/pfeToSed/calcToSed, SedDenitrif, OtherRemin, the
#####      bottom-cell coupling identity, and the assembled floor-cell O₂ consumption
#####   7. +cocco config — Eppley temp / VCO2 / picpoc + J_cocco*, freshly-built cocco baseline
#####
##### Grid: Oceananigans indexes k bottom→top, MARBL lev surface→bottom ⇒ oce k ↔ MARBL lev = nlev-k+1.
##### PAR is MARBL's PAR_avg, decomposed into ice-radiation sub-columns (FRACR_BIN/QSW_BIN) + the exact
##### interface PAR so the nonlinear Geider curve and the nitrification light-taper reproduce MARBL exactly.
#####
##### Out-of-scope terms belong to unimplemented phases (stated, not fudged): implicit ballast remin
##### (J_DOCr/DONr/DOPr, ballast flux into J_DIC/J_Alk/J_PO4/J_SiO3) → Phase 11 (here POC_REMIN_DIC is
##### supplied as the ballast input for the O₂ and sediment sections).
#####

using NCDatasets, Printf, Test, Oceananigans, OceanBioME
using Oceananigans.Fields: CenterField, Field, ConstantField, set!
using Oceananigans.Grids: Center, Face
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers, update_biogeochemical_state!
using Oceananigans.Units: day


include(joinpath(@__DIR__, "marbl_names.jl"))
include(joinpath(@__DIR__, "marbl_baselines.jl"))
const MP = @__MODULE__
require_baselines("compare_marbl") || exit(0)

include(joinpath(@__DIR__, "general_configs.jl"))   # the Phase 10 general N_A/N_Z configs
const ICM = OceanBioME.Models.NutrientsPlanktonDetritusModels.InorganicCarbonModels
using .ICM: biological_calcium_carbonate_precipitation, particulate_calcium_carbonate_production
using OceanBioME: NutrientsPlanktonDetritus, Nutrients, NitrateAmmonia, PO₄, Fe, Si, ExplicitCalciumCarbonate,
                  PrescribedPhotosyntheticallyActiveRadiation
using OceanBioME.Models.CarbonChemistryModel: CarbonChemistry, carbonate_concentration

# `carbon_dioxide_and_carbonate` is gone: CO₂aq is now a `CarbonChemistry` output option.
# Keep the old two-value form locally so the diagnostics below read unchanged.
carbon_dioxide_and_carbonate(cc; kwargs...) =
    (cc(; kwargs..., output = Val(:CO₂)), carbonate_concentration(cc; kwargs...))
const CCM = OceanBioME.Models.CarbonChemistryModel

# tiny median (Statistics is not a dep of the comparison env)
_median(v) = (s = sort(collect(v)); n = length(s); iseven(n) ? (s[n÷2] + s[n÷2+1]) / 2 : s[n÷2+1])

# ---- constant-ρ carbon chemistry + KF harness (test-only; do NOT touch OceanBioME defaults) ----
const CONST_RHO_CC = CarbonChemistry(; density_function = (args...) -> 1026.0)
# MARBL's KF: Dickson & Riley (1979), free scale → total via (1 + Sₜ/Kₛ); T in Kelvin.
struct DicksonRileyKF{PC}
    pressure_correction :: PC
end
@inline (c::DicksonRileyKF)(T, S, Is, KS; P = nothing) =
    c.pressure_correction(T, P) *
    exp(1590.2 / T - 12.641 + 1.525 * sqrt(Is) + log(1 - 0.001005 * S)
        + log(1 + (0.14 / 96.062 * S / 1.80655) / KS))
const CC_PEREZ = CarbonChemistry(; density_function = (args...) -> 1026.0,
    fluoride = CCM.KF(; constant = -9.68, inverse_T = 874.0, sqrt_S = 0.111))          # OceanBioME default
const CC_DR = CarbonChemistry(; density_function = (args...) -> 1026.0,
    fluoride = DicksonRileyKF(CCM.PressureCorrection(; a₀=-9.78, a₁=-0.0090, a₂=-0.000942, b₀=-0.00391, b₁=0.000054)))
# K1/K2 pressure correction with MARBL's gas-constant literal R=83.1451 (default OceanBioME 83.14472; both
# round the true 83.14463). Only demonstrates BIT-IDENTITY of the speciation when MARBL's pH is fed.
const CC_MARBL_R = CarbonChemistry(; density_function = (args...) -> 1026.0,
    carbonic_acid = (
        K1 = CCM.K1(; pressure_correction = CCM.PressureCorrection(; a₀=-25.50, a₁=0.1271, a₂=0.0, b₀=-0.00308, b₁=0.0000877, R=83.1451)),
        K2 = CCM.K2(; pressure_correction = CCM.PressureCorrection(; a₀=-15.82, a₁=-0.0219, a₂=0.0, b₀=0.00113, b₁=-0.0001475, R=83.1451))))

# ⭐⭐ THE MARBL-MATCHED CarbonChemistry — use THIS whenever testing carbonate against MARBL. An INDEPENDENT
# solve with it reproduces MARBL's H2CO3/CO3/pH to ~1e-13 (median) / ~1e-11 (max) — the Julia-vs-Fortran exp/log
# floor. ALL FOUR are required (drop one and it degrades — see project_carbonate_solver_comparison memory):
#   • DRTSafeSolver: MARBL's drtsafe (step-stop 1e-10) with MARBL's COLD bracket pH[6,9] ⇒ H_lower=1e-9, H_upper=1e-6.
#     ⚠️ the symmetric auto-bracket pH[7,9] MISSES deep pH<7 cells → 2.3e-4; only [6,9] reproduces MARBL.
#   • Dickson&Riley KF (default Perez&Fraga floors at 2.89e-7 — the KF was the whole gap).
#   • K1/K2 pressure R = 83.1451.   • sulfate 96.062 (source) + constant ρ = 1026.
const CC_MARBL = CarbonChemistry(; density_function = (args...) -> 1026.0,
    solver = OceanBioME.DRTSafeSolver(; max_iters = 100, bracket_grows = 3, atol = 1e-10, H_lower = 1e-9, H_upper = 1e-6),
    fluoride = DicksonRileyKF(CCM.PressureCorrection(; a₀=-9.78, a₁=-0.0090, a₂=-0.000942, b₀=-0.00391, b₁=0.000054)),
    carbonic_acid = (
        K1 = CCM.K1(; pressure_correction = CCM.PressureCorrection(; a₀=-25.50, a₁=0.1271, a₂=0.0, b₀=-0.00308, b₁=0.0000877, R=83.1451)),
        K2 = CCM.K2(; pressure_correction = CCM.PressureCorrection(; a₀=-15.82, a₁=-0.0219, a₂=0.0, b₀=0.00113, b₁=-0.0001475, R=83.1451))))

# hoisted for the array-reconstruction sections (a struct cannot be defined inside a let/function)
struct CheckRow; name::String; maxrel::Float64; n::Int; pass::Bool; dr::Float64; nd::Int; caught::Bool; end

# =====================================================================================================
# Shared preamble — base CESM2.1 baseline, one model on MARBL's exact column
# =====================================================================================================
const HIST = joinpath(MARBL_ROOT, "MARBL", "tests", "input_files", "baselines",
                      "call_compute_subroutines.history.nc")
const ICF  = joinpath(MARBL_ROOT, "MARBL", "tests", "input_files", "initial_conditions",
                      "call_compute_subroutines.20190718.nc")
const FILL = 9.96920996838687e36
ds = NCDataset(HIST)
r2(n) = permutedims(Array{Float64}(coalesce.(ds[n][:, :], FILL)))   # (col, lev)
fin(x) = (x != FILL) & isfinite(x)
zt = Array(ds["zt"][:]) / 100; zw = Array(ds["zw"][:]) / 100         # cm → m
ncols = size(r2("insitu_temp"), 1); nlev = length(zt)

zfaces = vcat(-reverse(zw), 0.0)
grid = RectilinearGrid(size = (ncols, 1, nlev), x = (0, ncols), y = (0, 1), z = zfaces,
                       topology = (Periodic, Periodic, Bounded))
flip(A) = (B = fill(0.0, ncols, 1, nlev);
           for c in 1:ncols, k in 1:nlev; v = A[c, nlev-k+1]; B[c,1,k] = fin(v) ? v : 0.0; end; B)

# ---- salinity + MARBL's exact prescribed pressure + ice-radiation bins (IC file) ----
icds = NCDataset(ICF)
Sraw    = permutedims(Array{Float64}(coalesce.(icds["salinity"][:, :], FILL)))
Praw_ic = permutedims(Array{Float64}(coalesce.(icds["pressure"][:, :], FILL)))   # bars, MARBL's exact field
FRACR = Array{Float64}(icds["FRACR_BIN"][:, :]); QSW = Array{Float64}(icds["QSW_BIN"][:, :])   # (nbin, col)
close(icds)
Smed = _median(filter(fin, Sraw))
salinity = (10 < Smed < 45) ? Sraw : Sraw .* 1000.0                              # apply scale_factor if packed

# ---- build the base current-config model (ExplicitCalciumCarbonate + MARBLOxygen), PAR prescribed ----
PARfield = CenterField(grid); light = PrescribedPhotosyntheticallyActiveRadiation(PARfield)
Sfield   = CenterField(grid)
bgc0 = NutrientsPlanktonDetritus(grid;
        nutrients = Nutrients(NitrateAmmonia(; nitrification_rate = 0.0), PO₄, Fe, Si),
        plankton  = MARBLPlankton(), detritus = MARBLDetritus(grid; sinking_speed = 10.0),
        inorganic_carbon = ExplicitCalciumCarbonate(grid; calcium_carbonate_dissolution_rate = 0.03/day, carbon_chemistry = CONST_RHO_CC),
        oxygen = MARBLOxygen(), light_attenuation = light)
model = NonhydrostaticModel(grid; tracers = PHYSICS_TRACERS, biogeochemistry = bgc0, buoyancy = nothing,
                            auxiliary_fields = (S = Sfield,))
bgc = model.biogeochemistry.underlying_biogeochemistry
clock = model.clock

# Constructing a MARBLPlankton is what DEFINES the tracer tendency methods for its PFT set (the names are
# only known then — see manifest_marbl_plankton!). A model built *inside* a function therefore calls them at
# a stale world age, where they are invisible and the abstraction's zero fallback (NPD L20) silently wins.
# The base config above is built at top level, so it is fine; `section_cocco` and `section_general` build
# their models inside the function, so manifest those PFT sets here (memoised — the later constructions
# reuse them).
const DUMMY_GRID = RectilinearGrid(size = (1, 1, 1), extent = (1, 1, 1))
MARBLPlankton_cocco(DUMMY_GRID)
general2_plankton()
general5_plankton(DUMMY_GRID)

marbl_of = Dict(:NO₃=>"NO3", :NH₄=>"NH4", :PO₄=>"PO4", :Fe=>"Fe", :Si=>"SiO3", :DIC=>"DIC", :Alk=>"ALK",
                :O₂=>"O2", :T=>"insitu_temp")
zerofill = Set((:POP, :PFe, :bSi, :POC, :CaCO₃))
tracers = model.tracers
for name in (required_biogeochemical_tracers(bgc)..., PHYSICS_TRACERS...)
    name in zerofill && continue
    set!(tracers[name], flip(r2(get(marbl_of, name, ascii_name(name)))))
end
set!(PARfield, flip(r2("PAR_avg")))
set!(Sfield,   flip(salinity))

# ---- ice-radiation sub-columns (φ, r) + exact mean-column interface PAR ----
nbin = size(FRACR, 1)
subcols_of(c) = (φ = @view FRACR[:, c]; meanQSW = sum(φ .* @view QSW[:, c]);
                 meanQSW > 0 ? ntuple(j -> (φ[j], (@view(QSW[:, c]) ./ meanQSW)[j]), nbin) :
                               ntuple(j -> (j == 1 ? 1.0 : 0.0, 0.0), nbin))
Δz_m = [zw_bot - zw_top for (zw_top, zw_bot) in zip(vcat(0.0, zw)[1:nlev], zw)]
KPARdz(chl, dz_m) = (w = max(chl, 0.02); (w < 0.13224 ? 0.000919 * w^0.3536 : 0.001131 * w^0.4562) * (dz_m * 100))
function build_PARiface(totalChl, PARavg, gridx)
    Pif = Field{Center, Center, Face}(gridx)
    for c in 1:ncols
        Iface = fill(0.0, nlev + 1)
        if fin(PARavg[c,1]) && PARavg[c,1] > 0
            K1 = KPARdz(totalChl[c,1], Δz_m[1])
            Iface[1] = PARavg[c,1] * K1 / (1 - exp(-K1))
            for l in 1:nlev; Iface[l+1] = Iface[l] * exp(-KPARdz(totalChl[c,l], Δz_m[l])); end
        end
        for kf in 1:(nlev + 1); @inbounds Pif[c, 1, kf] = Iface[nlev - kf + 2]; end
    end
    Pif
end
PARiface = build_PARiface(r2("spChl") .+ r2("diatChl") .+ r2("diazChl"), r2("PAR_avg"), grid)
auxs = ntuple(c -> (PAR = SubcolumnPAR(PARfield, map(first, subcols_of(c)), map(last, subcols_of(c))),
                    PAR_interface = PARiface), ncols)

# ---- masks (present-PFT / ghost-grazing exclusion) ----
Tm = r2("insitu_temp")
active(c, lev) = fin(Tm[c, lev])
bottomlev(c) = maximum(lev for lev in 1:nlev if active(c, lev))   # exclude the sediment-flux cell
pbase = MARBLPlankton()
A  = (sp = pbase.autotrophs.sp, diat = pbase.autotrophs.diat, diaz = pbase.autotrophs.diaz)
Cf = Dict(:sp => r2("spC"), :diat => r2("diatC"), :diaz => r2("diazC"))
# NB: `domcleanmask` (which used to drop cells where MARBL's `autotroph_zero_consistency_enforce` had
# silently killed a PFT) is GONE — the model now reproduces that zeroing internally (masked PFT-tracer
# reads, MARBLPlankton.jl), so those cells MATCH and every real-method comparison runs over the whole
# column. Same reason `section_general`'s input-state enforcement is no longer needed.

# ---- model-driven comparison engine (oce k ↔ marbl lev = nlev-k+1) ----
function mcompare(pred, truth; colmask = c -> true, cellmask = (c,lev) -> true, floor_frac = 1e-9)
    keep(c, lev) = colmask(c) && cellmask(c, lev) && active(c, lev) && lev != bottomlev(c) && fin(truth[c,lev])
    scale = 0.0
    for c in 1:ncols, lev in 1:nlev; keep(c,lev) && (scale = max(scale, abs(truth[c,lev]))); end
    thr = floor_frac * scale
    mr = 0.0; n = 0; worst = (0,0,0.0,0.0)
    for c in 1:ncols, k in 1:nlev
        lev = nlev - k + 1
        keep(c, lev) || continue
        abs(truth[c,lev]) ≤ thr && continue
        p = pred(c, k, auxs[c]); isfinite(p) || continue
        n += 1
        r = abs(p - truth[c,lev]) / (abs(truth[c,lev]) + 1e-14)
        r > mr && (mr = r; worst = (c, lev, p, truth[c,lev]))
    end
    (mr, n, worst)
end
mreport(name, mr, n, worst; tol) = begin
    pass = mr ≤ tol
    tag = mr ≤ 1e-10 ? "MACHINE-PRECISION ✅" : pass ? "PASS ✅" : "FAIL ❌"
    @printf("  %-18s max-rel=%.2e over %4d cells  %s\n", name, mr, n, tag)
    pass || @printf("      worst @ (c=%d,lev=%d) ours=% .6e marbl=% .6e\n", worst...)
    pass
end

# ---- shared array-reconstruction comparison (records dynamic range + perturbation-catchability) ----
function acompare!(results, name, predict, truth; atol = 1e-14, mask = trues(size(truth)))
    mr = 0.0; n = 0; vals = Float64[]; allclose = true; caught = false; worst = (0,0.0,0.0)
    for I in eachindex(truth)
        (mask[I] && fin(truth[I])) || continue
        pr = predict[I]; isfinite(pr) || continue
        n += 1; r = abs(pr - truth[I]) / (abs(truth[I]) + atol)
        r > mr && (mr = r; worst = (I, pr, truth[I]))
        # `atol` floors BOTH the reported metric (above) and the pass/fail — otherwise two values that are
        # both numerical dust (e.g. ±1e-26 for a field whose real scale is 1e-9) are ranked machine-precise
        # by `r` yet declared "not close" by a bare isapprox, and the row prints a tiny max-rel but FAILs.
        allclose &= isapprox(pr, truth[I]; atol = atol); caught |= !isapprox(pr*(1+1e-3), truth[I]; atol = atol)
        truth[I] != 0 && push!(vals, abs(truth[I]))
    end
    dr = isempty(vals) ? 0.0 : maximum(vals)/minimum(vals)
    nd = length(unique(round.(vals, sigdigits=6)))
    push!(results, CheckRow(name, mr, n, allclose, dr, nd, caught))
    allclose || (@printf("    [worst %-26s] idx=%s pred=% .6e marbl=% .6e\n", name,
                         Tuple(CartesianIndices(truth)[worst[1]]), worst[2], worst[3]))
end
function report_rows(results)
    pass = true
    for r in results
        note = !r.caught ? "⚠ not-catchable" : (r.dr < 2 || r.nd < 5) ? "weak (degenerate)" : "strong"
        @printf("  %-32s max-rel=%.2e over %4d cells  %s  %s\n", r.name, r.maxrel, r.n, r.pass ? "PASS ✅" : "FAIL ❌", note)
        p = r.pass; @test p; pass &= p
    end
    pass
end

# =====================================================================================================
# Section 1 — per-RATE terms + grazing network (phase-3/4 milestone) + photoC sub-column reconstruction
# =====================================================================================================
function section_rates()
    println("\n═══ 1. per-rate terms + grazing network vs MARBL ═══")
    T = r2("insitu_temp"); PAR_avg = r2("PAR_avg")
    NO3=r2("NO3"); NH4=r2("NH4"); PO4=r2("PO4"); FeA=r2("Fe"); SiO3=r2("SiO3"); DOP=r2("DOP")
    spC=r2("spC"); spP=r2("spP"); spFe=r2("spFe"); spCaCO3=r2("spCaCO3"); spChl=r2("spChl")
    diatC=r2("diatC"); diatP=r2("diatP"); diatFe=r2("diatFe"); diatSi=r2("diatSi")
    diazC=r2("diazC"); diazP=r2("diazP"); diazFe=r2("diazFe"); zooC=r2("zooC")
    M = Dict(n => r2(n) for n in
        ("photoC_sp","photoC_diat","photoC_diaz","photoFe_sp","photoFe_diat","photoFe_diaz",
         "photoNO3_sp","photoNO3_diat","photoNO3_diaz","photoNH4_sp","photoNH4_diat","photoNH4_diaz",
         "PO4_sp_uptake","PO4_diat_uptake","PO4_diaz_uptake","DOP_sp_uptake","DOP_diat_uptake","DOP_diaz_uptake",
         "sp_Qp","diat_Qp","diaz_Qp","sp_loss","diat_loss","diaz_loss",
         "sp_loss_poc","sp_loss_doc","diat_loss_doc","diaz_loss_doc","sp_agg","diat_agg","diaz_agg",
         "sp_CaCO3_form","diat_bSi_form","diaz_Nfix",
         "graze_sp","graze_diat","graze_diaz","graze_sp_zoo","graze_diat_zoo","graze_diaz_zoo",
         "graze_sp_poc","graze_diat_poc","graze_diaz_poc","graze_sp_doc","graze_diat_doc","graze_diaz_doc",
         "zoo_loss","zoo_loss_poc","zoo_loss_doc","x_graze_zoo","graze_zoo","J_zooC"))
    p = MARBLPlankton()
    Aa = (sp = p.autotrophs.sp, diat = p.autotrophs.diat, diaz = p.autotrophs.diaz)
    Q = p.nitrogen_to_carbon; εC = p.concentration_regularisation
    rFe = p.iron_quota_threshold_factor; rSi = p.silicon_quota_threshold_factor
    khs(s) = Aa[s].nutrient_half_saturations
    depth_m(lev) = abs(zt[lev])                     # zt already in metres (preamble /100)
    function Pprime(s, C, Tloc, lev)
        a = Aa[s]; z2, z1 = p.loss_threshold_zero_depth, p.loss_threshold_full_depth
        f_thres = clamp((z2 - depth_m(lev)) / (z2 - z1), 0.0, 1.0)
        thres = Tloc < a.temperature_threshold ? a.cold_loss_threshold_concentration : a.loss_threshold_concentration
        max(C - f_thres * thres, 0.0)
    end
    function fnut_marbl(s, no3, nh4, po4, dop, fe, sio3)
        a = Aa[s]; k = khs(s); x = no3/k.nitrate; y = nh4/k.ammonia
        VN = s === :diaz ? 1.0 : (x + y)/(1 + x + y); VFe = fe/(fe + k.iron)
        u = po4/k.phosphate; v = dop/k.DOP; VP = (u + v)/(1 + u + v)
        f = min(VN, VFe, VP); s === :diat ? min(f, sio3/(sio3 + k.silicate)) : f
    end
    Csp(s) = s === :sp ? spC : (s === :diat ? diatC : diazC)
    act = fin.(T)
    pred = Dict{String,Matrix{Float64}}(); allocp(name) = (pred[name] = fill(NaN, ncols, nlev))
    for arr in ("sp_Qp","diat_Qp","diaz_Qp","sp_loss","diat_loss","diaz_loss","sp_loss_poc","sp_loss_doc",
                "diat_loss_doc","diaz_loss_doc","sp_agg","diat_agg","diaz_agg","gQfe_sp","gQfe_diat","gQfe_diaz",
                "gQp_sp","gQp_diat","gQp_diaz","gQsi_diat","Nuptake_sp","Nuptake_diat","Nsplit_sp","Nsplit_diat",
                "Nsplit_diaz","Nfix_diaz","CaCO3_sp","fnut_sp_caco3"); allocp(arr); end
    for c in 1:ncols, k in 1:nlev
        Tl = T[c,k]
        pred["sp_Qp"][c,k]=spP[c,k]/(spC[c,k]+εC); pred["diat_Qp"][c,k]=diatP[c,k]/(diatC[c,k]+εC); pred["diaz_Qp"][c,k]=diazP[c,k]/(diazC[c,k]+εC)
        for s in (:sp,:diat,:diaz)
            a = Aa[s]; C = Csp(s)[c,k]; Pp = Pprime(s, C, Tl, k)
            Tf = temperature_scaling(a.temperature_sensitivity, a.temperature_reference, Tl)
            pred["$(s)_loss"][c,k] = a.linear_mortality_rate * Pp * Tf
            pred["$(s)_agg"][c,k]  = aggregation_rate(a.quadratic_mortality_rate, p.aggregation_exponent,
                                                      a.minimum_aggregation_rate, a.maximum_aggregation_rate, Pp)
        end
        QCa = min(spCaCO3[c,k] / (spC[c,k] + εC), p.maximum_calcite_quota)
        pred["sp_loss_poc"][c,k]   = QCa * pred["sp_loss"][c,k]
        pred["sp_loss_doc"][c,k]   = (1 - p.labile_fraction) * (pred["sp_loss"][c,k] - pred["sp_loss_poc"][c,k])
        pred["diat_loss_doc"][c,k] = (1 - p.labile_fraction) * pred["diat_loss"][c,k]
        pred["diaz_loss_doc"][c,k] = (1 - p.labile_fraction) * pred["diaz_loss"][c,k]
        for s in (:sp,:diat,:diaz)
            a = Aa[s]; pc = M["photoC_$s"][c,k]
            gQfe = growth_iron_quota(a.growth_iron_quota_reference, a.minimum_iron_quota, rFe, khs(s).iron, FeA[c,k])
            gQp  = growth_phosphorus_quota(p.phosphate_quota_slope, p.phosphate_quota_intercept, p.phosphate_to_carbon_maximum, PO4[c,k])
            pred["gQfe_$s"][c,k] = pc * gQfe; pred["gQp_$s"][c,k] = pc * gQp
        end
        pred["gQsi_diat"][c,k] = M["photoC_diat"][c,k] *
            growth_silicon_quota(p.silicon_quota_reference, p.maximum_silicon_quota, p.minimum_silicon_quota,
                                 rFe, rSi, khs(:diat).iron, khs(:diat).silicate, FeA[c,k], SiO3[c,k])
        pred["Nuptake_sp"][c,k] = M["photoC_sp"][c,k]*Q; pred["Nuptake_diat"][c,k] = M["photoC_diat"][c,k]*Q
        for s in (:sp,:diat,:diaz)
            k_ = khs(s); xr = NO3[c,k]/k_.nitrate; yr = NH4[c,k]/k_.ammonia; pred["Nsplit_$s"][c,k] = xr/(xr+yr)
        end
        pred["Nfix_diaz"][c,k] = nitrogen_fixation(p.nitrogen_fixation_ratio, Q, M["photoC_diaz"][c,k], M["photoNO3_diaz"][c,k], M["photoNH4_diaz"][c,k])
        fn = fnut_marbl(:sp, NO3[c,k], NH4[c,k], PO4[c,k], DOP[c,k], FeA[c,k], SiO3[c,k])
        pred["CaCO3_sp"][c,k] = calcite_formation(p.calcite_production_fraction, p.maximum_calcifying_fraction,
                                    p.calcite_bloom_threshold, p.calcite_temperature_threshold, p.calcite_temperature_minimum,
                                    M["photoC_sp"][c,k], fn, Tl, spC[c,k])
        (Tl >= p.calcite_temperature_threshold) && (spC[c,k] < p.calcite_bloom_threshold) && (M["photoC_sp"][c,k] > 0) && (pred["fnut_sp_caco3"][c,k] = fn)
    end
    pres(s) = act .& (Csp(s) .> Aa[s].loss_threshold_concentration)
    hasC(s) = act .& (M["photoC_$s"] .> 0)
    calc_nomod = act .& (T .>= p.calcite_temperature_threshold) .& (spC .< p.calcite_bloom_threshold) .& (M["photoC_sp"] .> 0)
    fnut_truth = fill(FILL, ncols, nlev)
    for c in 1:ncols, k in 1:nlev
        calc_nomod[c,k] && (fnut_truth[c,k] = sqrt(M["sp_CaCO3_form"][c,k] / (p.calcite_production_fraction * M["photoC_sp"][c,k])))
    end
    Nsplit_truth(s) = M["photoNO3_$s"] ./ (M["photoNO3_$s"] .+ M["photoNH4_$s"])
    Nsplit_mask(s)  = act .& ((M["photoNO3_$s"] .+ M["photoNH4_$s"]) .> 1e-30)

    results = CheckRow[]; C! = (a...; kw...) -> acompare!(results, a...; kw...)
    C!("sp_Qp  (standing P:C)",  pred["sp_Qp"],   M["sp_Qp"];   mask = pres(:sp))
    C!("diat_Qp (standing P:C)", pred["diat_Qp"], M["diat_Qp"]; mask = pres(:diat))
    C!("diaz_Qp (standing P:C)", pred["diaz_Qp"], M["diaz_Qp"]; mask = pres(:diaz))
    C!("sp_loss (mort·P′·Tf)",   pred["sp_loss"],   M["sp_loss"];   mask = pres(:sp))
    C!("diat_loss (mort·P′·Tf)", pred["diat_loss"], M["diat_loss"]; mask = pres(:diat))
    C!("diaz_loss (mort·P′·Tf)", pred["diaz_loss"], M["diaz_loss"]; mask = pres(:diaz))
    C!("sp_agg (aggregation)",   pred["sp_agg"],   M["sp_agg"];   mask = pres(:sp))
    C!("diat_agg (aggregation)", pred["diat_agg"], M["diat_agg"]; mask = pres(:diat))
    C!("diaz_agg (aggregation)", pred["diaz_agg"], M["diaz_agg"]; mask = pres(:diaz))
    C!("sp loss→POC (ballast)",  pred["sp_loss_poc"],   M["sp_loss_poc"];   mask = pres(:sp))
    C!("sp loss→DOC (labile)",   pred["sp_loss_doc"],   M["sp_loss_doc"];   mask = pres(:sp))
    C!("diat loss→DOC",          pred["diat_loss_doc"], M["diat_loss_doc"]; mask = pres(:diat))
    C!("diaz loss→DOC",          pred["diaz_loss_doc"], M["diaz_loss_doc"]; mask = pres(:diaz))
    C!("P uptake sp = pC·gQp",   pred["gQp_sp"],   M["PO4_sp_uptake"]  .+ M["DOP_sp_uptake"];   mask = hasC(:sp))
    C!("P uptake diat = pC·gQp", pred["gQp_diat"], M["PO4_diat_uptake"].+ M["DOP_diat_uptake"]; mask = hasC(:diat))
    C!("P uptake diaz = pC·gQp", pred["gQp_diaz"], M["PO4_diaz_uptake"].+ M["DOP_diaz_uptake"]; mask = hasC(:diaz))
    C!("photoFe sp = pC·gQfe",   pred["gQfe_sp"],   M["photoFe_sp"];   mask = hasC(:sp))
    C!("photoFe diat = pC·gQfe", pred["gQfe_diat"], M["photoFe_diat"]; mask = hasC(:diat))
    C!("photoFe diaz = pC·gQfe", pred["gQfe_diaz"], M["photoFe_diaz"]; mask = hasC(:diaz))
    C!("photoSi diat = pC·gQsi", pred["gQsi_diat"], M["diat_bSi_form"]; mask = hasC(:diat))
    C!("N uptake sp = pC·Q",     pred["Nuptake_sp"],  M["photoNO3_sp"]  .+ M["photoNH4_sp"];   mask = hasC(:sp))
    C!("N uptake diat = pC·Q",   pred["Nuptake_diat"],M["photoNO3_diat"].+ M["photoNH4_diat"]; mask = hasC(:diat))
    C!("N split sp (Monod)",     pred["Nsplit_sp"],  Nsplit_truth(:sp);   mask = Nsplit_mask(:sp))
    C!("N split diat (Monod)",   pred["Nsplit_diat"],Nsplit_truth(:diat); mask = Nsplit_mask(:diat))
    C!("N split diaz (Monod)",   pred["Nsplit_diaz"],Nsplit_truth(:diaz); mask = Nsplit_mask(:diaz))
    C!("diaz_Nfix (§10)",        pred["Nfix_diaz"],  M["diaz_Nfix"];      mask = hasC(:diaz))
    C!("sp_CaCO3_form (§9.1)",   pred["CaCO3_sp"],   M["sp_CaCO3_form"];  mask = hasC(:sp))
    C!("f_nut(sp) via CaCO3",    pred["fnut_sp_caco3"], fnut_truth;       mask = calc_nomod)

    # --- Phase-4 grazing network (§13, §14) ---
    zp = p.zooplankton.zoo
    function Zprime_zoo(zc, lev)
        z2, z1 = p.zoo_loss_threshold_zero_depth, p.zoo_loss_threshold_full_depth
        f_thres = clamp((z2 - depth_m(lev)) / (z2 - z1), 0.0, 1.0); max(zc - f_thres * zp.loss_threshold_concentration, 0.0)
    end
    for arr in ("graze_sp","graze_diat","graze_diaz","graze_sp_zoo","graze_diat_zoo","graze_diaz_zoo",
                "graze_sp_poc","graze_diat_poc","graze_diaz_poc","graze_sp_doc","graze_diat_doc","graze_diaz_doc",
                "zoo_loss_p","zoo_loss_poc_p","zoo_loss_doc_p","x_graze_zoo_p","J_zooC_p"); allocp(arr); end
    εTinv = p.growth_rate_regularisation
    for c in 1:ncols, k in 1:nlev
        Tl = T[c,k]; zc = zooC[c,k]; Tf_z = temperature_scaling(zp.temperature_sensitivity, zp.temperature_reference, Tl)
        ag = Dict{Symbol,Float64}()
        for s in (:sp,:diat,:diaz)
            g = zp.grazing[s]; Pp = Pprime(s, Csp(s)[c,k], Tl, k)
            graze = grazing_rate(g.maximum_grazing_rate, g.grazing_half_saturation, g.sigmoidal, Tf_z, zc, Pp)
            ag[s] = graze; pred["graze_$s"][c,k] = graze
            pred["graze_$(s)_zoo"][c,k] = g.fraction_to_zooplankton * graze
            pred["graze_$(s)_doc"][c,k] = g.fraction_to_dissolved * graze
            if s === :sp
                QCa = min(spCaCO3[c,k] / (spC[c,k] + εC), p.maximum_calcite_quota)
                frac = calcifier_graze_poc_fraction(p.calcite_ballast_minimum, QCa, p.small_phyto_poc_factor, Pp, p.grazed_small_phyto_poc_limit)
                pred["graze_sp_poc"][c,k] = frac * graze
            else
                pred["graze_$(s)_poc"][c,k] = g.fraction_to_particulate * graze
            end
        end
        Zp = Zprime_zoo(zc, k); zl = zoo_loss_rate(zp.quadratic_mortality_rate, zp.linear_mortality_rate, p.zoo_aggregation_exponent, Zp, Tf_z)
        ε = εC * εTinv; wsp, wdiat, wdiaz = ag[:sp] + ε, ag[:diat] + ε, ag[:diaz] + ε
        fzd = (zp.grazing.sp.detritus_fraction*wsp + zp.grazing.diat.detritus_fraction*wdiat + zp.grazing.diaz.detritus_fraction*wdiaz) / (wsp + wdiat + wdiaz)
        pred["zoo_loss_p"][c,k] = zl; pred["zoo_loss_poc_p"][c,k] = fzd * zl
        pred["zoo_loss_doc_p"][c,k] = (1 - p.labile_fraction) * (1 - fzd) * zl
        xg = pred["graze_sp_zoo"][c,k] + pred["graze_diat_zoo"][c,k] + pred["graze_diaz_zoo"][c,k]
        pred["x_graze_zoo_p"][c,k] = xg; pred["J_zooC_p"][c,k] = xg - zl
    end
    grpres(s) = act .& (M["graze_$s"] .> 0); zlpres = act .& (M["zoo_loss"] .> 0)
    graze_active = act .& (pres(:sp) .| pres(:diat) .| pres(:diaz))
    C!("graze sp (auto_graze)",  pred["graze_sp"],   M["graze_sp"];   mask = grpres(:sp))
    C!("graze diat (auto_graze)",pred["graze_diat"], M["graze_diat"]; mask = grpres(:diat))
    C!("graze diaz (auto_graze)",pred["graze_diaz"], M["graze_diaz"]; mask = grpres(:diaz))
    C!("graze sp →zoo",   pred["graze_sp_zoo"],   M["graze_sp_zoo"];   mask = grpres(:sp))
    C!("graze diat→zoo",  pred["graze_diat_zoo"], M["graze_diat_zoo"]; mask = grpres(:diat))
    C!("graze diaz→zoo",  pred["graze_diaz_zoo"], M["graze_diaz_zoo"]; mask = grpres(:diaz))
    C!("graze sp →POC (ballast)", pred["graze_sp_poc"],   M["graze_sp_poc"];   mask = grpres(:sp))
    C!("graze diat→POC",          pred["graze_diat_poc"], M["graze_diat_poc"]; mask = grpres(:diat))
    C!("graze diaz→POC",          pred["graze_diaz_poc"], M["graze_diaz_poc"]; mask = grpres(:diaz))
    C!("graze sp →DOC",   pred["graze_sp_doc"],   M["graze_sp_doc"];   mask = grpres(:sp))
    C!("graze diat→DOC",  pred["graze_diat_doc"], M["graze_diat_doc"]; mask = grpres(:diat))
    C!("graze diaz→DOC",  pred["graze_diaz_doc"], M["graze_diaz_doc"]; mask = grpres(:diaz))
    C!("zoo_loss (mort₂Z′^e+mortZ′)Tf", pred["zoo_loss_p"],     M["zoo_loss"];     mask = zlpres)
    C!("zoo_loss→POC (f_zoo_detr)",     pred["zoo_loss_poc_p"], M["zoo_loss_poc"]; mask = zlpres .& graze_active)
    C!("zoo_loss→DOC (labile)",         pred["zoo_loss_doc_p"], M["zoo_loss_doc"]; mask = zlpres .& graze_active)
    C!("x_graze_zoo (Σ →zoo)",          pred["x_graze_zoo_p"],  M["x_graze_zoo"];  mask = graze_active)
    C!("J_zooC = x_graze − zoo_loss",   pred["J_zooC_p"],       M["J_zooC"];       mask = graze_active)
    pass = report_rows(results)

    # --- photoC vs MARBL (6 sub-column reconstruction) + negative control ---
    biomass(s) = s === :sp ? (spC, spChl) : s === :diat ? (diatC, r2("diatChl")) : (diazC, r2("diazChl"))
    function photoC_subcol(s, α)
        truth = M["photoC_$s"]; a = Aa[s]; C, Chl = biomass(s); mr = 0.0; n = 0
        for c in 1:ncols
            φ = @view FRACR[:, c]; meanQSW = sum(φ .* @view QSW[:, c]); r = meanQSW > 0 ? (@view(QSW[:, c]) ./ meanQSW) : zero(φ)
            for k in 1:nlev
                i = CartesianIndex(c, k); (act[i] && truth[i] > 0) || continue; n += 1
                Tf = temperature_scaling(a.temperature_sensitivity, a.temperature_reference, T[i])
                fn = fnut_marbl(s, NO3[i], NH4[i], PO4[i], DOP[i], FeA[i], SiO3[i])
                PCmax = max_specific_rate(a.photosynthesis_rate_reference, fn, Tf, T[i], a.temperature_threshold)
                coef = α * (Chl[i] / (C[i] + εC)) / (PCmax + p.growth_rate_regularisation); PCphoto = 0.0
                for j in eachindex(φ); PARj = PAR_avg[i] * r[j]; PARj > 0 && (PCphoto += φ[j] * PCmax * (1 - exp(-coef * PARj))); end
                mr = max(mr, abs(PCphoto * C[i] - truth[i]) / truth[i])
            end
        end
        (mr, n)
    end
    println("── photoC vs MARBL (reconstruct MARBL's 6 surface-light sub-columns) ──")
    for s in (:sp, :diat, :diaz)
        mr, n = photoC_subcol(s, Aa[s].initial_PI_slope)
        @printf("  photoC %-4s max-rel=%.2e over %4d cells  %s\n", s, mr, n, mr ≤ 1e-10 ? "MACHINE-PRECISION ✅" : "FAIL ❌")
        pp = mr ≤ 1e-10; @test pp; pass &= pp
    end
    bad, = photoC_subcol(:diat, Aa[:diat].initial_PI_slope * (0.39/0.28))     # negative control (the alphaPI bug)
    @printf("  [neg ctrl] diat alphaPI=0.39 → max-rel=%.2e caught? %s\n", bad, bad > 1e-3 ? "YES ✅" : "NO ❌")
    @test bad > 1e-3
    pass
end

# =====================================================================================================
# Section 2 — REAL METHOD  bgc(Val:X)  vs MARBL J_<tracer>  (autotroph quotas + DOM; sub-columns exact)
# =====================================================================================================
function section_realmethod()
    println("\n═══ 2. REAL METHOD  bgc(i,j,k,Val:X)  vs MARBL J_<tracer> ═══")
    set!(tracers[:POC], flip(fill(0.0, ncols, nlev)))
    autos = [(:spC,:sp),(:diatC,:diat),(:diazC,:diaz),(:spChl,:sp),(:diatChl,:diat),(:diazChl,:diaz),
             (:spP,:sp),(:diatP,:diat),(:diazP,:diaz),(:spFe,:sp),(:diatFe,:diat),(:diazFe,:diaz),
             (:diatSi,:diat),(:spCaCO₃,:sp)]
    presmask(s) = (c,lev) -> Cf[s][c,lev] > A[s].loss_threshold_concentration
    pass = true
    println("── autotroph quotas (machine precision) ──")
    for (nm, s) in autos
        r = mcompare((c,k,aux) -> bgc(c,1,k, grid, Val(nm), clock, tracers, aux), r2("J_" * ascii_name(nm));
                     cellmask = presmask(s), floor_frac = 0.0)
        p = mreport("J_" * ascii_name(nm), r...; tol = 1e-10); @test p; pass &= p
    end
    println("── semilabile DOM (machine precision; whole column — zeroing enforced in-model) ──")
    for nm in (:DOC, :DON, :DOP)
        r = mcompare((c,k,aux) -> bgc(c,1,k, grid, Val(nm), clock, tracers, aux), r2("J_" * ascii_name(nm));
                     floor_frac = 0.0)
        p = mreport("J_" * ascii_name(nm), r...; tol = 1e-10); @test p; pass &= p
    end
    pass
end

# =====================================================================================================
# Section 3 — ASSEMBLED J_<tracer> reconstructed from MARBL rate diagnostics (validates photoacc synthesis)
# =====================================================================================================
function section_jtend()
    println("\n═══ 3. ASSEMBLED J_<tracer> (reconstructed from MARBL rate diagnostics) ═══")
    T = r2("insitu_temp"); PAR_avg = r2("PAR_avg")
    NO3=r2("NO3"); NH4=r2("NH4"); PO4=r2("PO4"); FeA=r2("Fe"); SiO3=r2("SiO3"); DOP=r2("DOP")
    state = Dict(n => r2(n) for n in
      ("spC","spChl","spP","spFe","spCaCO3","diatC","diatChl","diatP","diatFe","diatSi","diazC","diazChl","diazP","diazFe"))
    M = Dict(n => r2(n) for n in
      ("photoC_sp","photoC_diat","photoC_diaz","photoFe_sp","photoFe_diat","photoFe_diaz",
       "PO4_sp_uptake","PO4_diat_uptake","PO4_diaz_uptake","DOP_sp_uptake","DOP_diat_uptake","DOP_diaz_uptake",
       "graze_sp","graze_diat","graze_diaz","sp_loss","diat_loss","diaz_loss","sp_agg","diat_agg","diaz_agg",
       "sp_CaCO3_form","diat_bSi_form",
       "DOC_prod","DON_prod","DOP_prod","DOC_remin","DON_remin","DOP_remin","DOP_loss_P_bal",
       "J_spC","J_diatC","J_diazC","J_spP","J_diatP","J_diazP","J_spChl","J_diatChl","J_diazChl",
       "J_spFe","J_diatFe","J_diazFe","J_diatSi","J_spCaCO3","J_DOC","J_DON","J_DOP"))
    Q  = pbase.nitrogen_to_carbon; εC = pbase.concentration_regularisation
    khs(s) = A[s].nutrient_half_saturations
    function fnut(s, no3, nh4, po4, dop, fe, sio3)
        a = A[s]; k = khs(s); x = no3/k.nitrate; y = nh4/k.ammonia
        VN = s === :diaz ? 1.0 : (x + y)/(1 + x + y); VFe = fe/(fe + k.iron)
        u = po4/k.phosphate; v = dop/k.DOP; VP = (u + v)/(1 + u + v)
        f = min(VN, VFe, VP); s === :diat ? min(f, sio3/(sio3 + k.silicate)) : f
    end
    Cf2(s)  = state[s === :sp ? "spC"   : s === :diat ? "diatC"   : "diazC"]
    Chlf(s) = state[s === :sp ? "spChl" : s === :diat ? "diatChl" : "diazChl"]
    function photoacc_recon(s, i, c)
        a = A[s]; C = Cf2(s)[i]; Chl = Chlf(s)[i]; θC = Chl / (C + εC)
        Tf = temperature_scaling(a.temperature_sensitivity, a.temperature_reference, T[i])
        fn = fnut(s, NO3[i], NH4[i], PO4[i], DOP[i], FeA[i], SiO3[i])
        PCmax = max_specific_rate(a.photosynthesis_rate_reference, fn, Tf, T[i], a.temperature_threshold)
        φ = @view FRACR[:, c]; meanQSW = sum(φ .* @view QSW[:, c]); r = meanQSW > 0 ? (@view(QSW[:, c]) ./ meanQSW) : zero(φ)
        photoacc = 0.0
        for j in eachindex(φ)
            PARj = PAR_avg[i] * r[j]
            if PARj > 0
                photoacc += φ[j] * photoacclimation_rate(a, pbase, PCmax, θC, Chl, PARj)
            end
        end
        photoacc
    end
    J = Dict{String,Matrix{Float64}}(); for name in keys(M); startswith(name, "J_") && (J[name] = fill(NaN, ncols, nlev)); end
    act = fin.(T)
    for c in 1:ncols, k in 1:nlev
        i = CartesianIndex(c, k); act[i] || continue
        for s in (:sp, :diat, :diaz)
            Cs = Cf2(s)[i]; auto_sum = M["graze_$s"][i] + M["$(s)_loss"][i] + M["$(s)_agg"][i]; fin(auto_sum) || continue
            invC = 1 / (Cs + εC)
            J["J_$(s)C"][i]   = M["photoC_$s"][i] - auto_sum
            J["J_$(s)P"][i]   = M["PO4_$(s)_uptake"][i] + M["DOP_$(s)_uptake"][i] - state["$(s)P"][i]*invC * auto_sum
            J["J_$(s)Fe"][i]  = M["photoFe_$s"][i] - state["$(s)Fe"][i]*invC * auto_sum
            J["J_$(s)Chl"][i] = photoacc_recon(s, i, c) - Chlf(s)[i]*invC * auto_sum
        end
        asd = M["graze_diat"][i] + M["diat_loss"][i] + M["diat_agg"][i]
        fin(asd) && (J["J_diatSi"][i] = M["diat_bSi_form"][i] - state["diatSi"][i]/(state["diatC"][i]+εC) * asd)
        ass = M["graze_sp"][i] + M["sp_loss"][i] + M["sp_agg"][i]
        fin(ass) && (J["J_spCaCO3"][i] = M["sp_CaCO3_form"][i] - state["spCaCO3"][i]/(state["spC"][i]+εC) * ass)
    end
    rC, rN, rP = 0.01, 0.0115, 0.003
    for c in 1:ncols, k in 1:nlev
        i = CartesianIndex(c, k); act[i] || continue
        J["J_DOC"][i] = (1-rC)*M["DOC_prod"][i] - M["DOC_remin"][i]
        J["J_DON"][i] = (1-rN)*M["DON_prod"][i] - M["DON_remin"][i]
        J["J_DOP"][i] = (1-rP)*M["DOP_prod"][i] - M["DOP_remin"][i] -
                        M["DOP_sp_uptake"][i] - M["DOP_diat_uptake"][i] - M["DOP_diaz_uptake"][i] - M["DOP_loss_P_bal"][i]
    end
    pres(s) = act .& (Cf2(s) .> A[s].loss_threshold_concentration)
    function check(name, mask)
        pred = J[name]; truth = M[name]; mr = 0.0; n = 0; caught = false; allclose = true
        for i in eachindex(truth)
            (mask[i] && fin(truth[i]) && isfinite(pred[i])) || continue
            n += 1; r = abs(pred[i]-truth[i])/(abs(truth[i])+1e-14); mr = max(mr, r)
            allclose &= isapprox(pred[i], truth[i]); caught |= !isapprox(pred[i]*(1+1e-3), truth[i])
        end
        tag = allclose ? (mr ≤ 1e-10 ? "MACHINE-PRECISION ✅" : "PASS ✅") : "FAIL ❌"
        @printf("  %-12s max-rel=%.2e over %4d cells  %s%s\n", name, mr, n, tag, caught ? "" : "  ⚠ not-catchable")
        allclose
    end
    pass = true
    for s in (:sp,:diat,:diaz); p = check("J_$(s)C", pres(s));   @test p; pass &= p; end
    for s in (:sp,:diat,:diaz); p = check("J_$(s)P", pres(s));   @test p; pass &= p; end
    for s in (:sp,:diat,:diaz); p = check("J_$(s)Chl", pres(s)); @test p; pass &= p; end
    for s in (:sp,:diat,:diaz); p = check("J_$(s)Fe", pres(s));  @test p; pass &= p; end
    p = check("J_diatSi", pres(:diat)); @test p; pass &= p
    p = check("J_spCaCO3", pres(:sp));  @test p; pass &= p
    for nm in ("J_DOC","J_DON","J_DOP"); p = check(nm, act); @test p; pass &= p; end
    pass
end

# =====================================================================================================
# Section 4 — DOM / POM / PFe rate terms (§15.1/§15.2/§16.1) reconstructed from MARBL diagnostics
# =====================================================================================================
function section_dompom()
    println("\n═══ 4. DOM / POM / PFe rate terms vs MARBL ═══")
    T = r2("insitu_temp"); PAR_avg = r2("PAR_avg")
    DOC=r2("DOC"); DON=r2("DON"); DOP=r2("DOP"); DOCr=r2("DOCr"); DONr=r2("DONr"); DOPr=r2("DOPr")
    spC=r2("spC"); spP=r2("spP"); diatC=r2("diatC"); diatP=r2("diatP"); diazC=r2("diazC"); diazP=r2("diazP")
    spFe=r2("spFe"); diatFe=r2("diatFe"); diazFe=r2("diazFe")
    M = Dict(n => r2(n) for n in
        ("DOC_prod","DON_prod","DOP_prod","DOC_remin","DON_remin","DOP_remin","DOCr_remin","DONr_remin","DOPr_remin",
         "POC_REMIN_DIC","POC_REMIN_DOCr","POP_REMIN_PO4","POP_REMIN_DOPr","PON_REMIN_NH4","PON_REMIN_DONr",
         "graze_sp","graze_diat","graze_diaz","sp_loss","diat_loss","diaz_loss","sp_agg","diat_agg","diaz_agg",
         "graze_sp_poc","graze_diat_poc","graze_diaz_poc","sp_loss_poc","diat_loss_poc","diaz_loss_poc",
         "graze_sp_zoo","graze_diat_zoo","graze_diaz_zoo","sp_loss_doc","diat_loss_doc","diaz_loss_doc",
         "graze_sp_doc","graze_diat_doc","graze_diaz_doc",
         "zoo_loss_doc","zoo_loss_poc","graze_zoo_doc","graze_zoo_poc","P_iron_PROD","Fe_scavenge"))
    p = MARBLPlankton(); d = MARBLDetritus(RectilinearGrid(size=(1,1,1), extent=(1,1,1)))
    Q = p.nitrogen_to_carbon; f_toDON = p.nitrogen_to_DON_fraction; f_toDOP = p.phosphorus_to_DOP_fraction
    Qp_zoo = p.zooplankton_phosphate_to_carbon; εC = p.concentration_regularisation; Qfe_zoo = p.zooplankton_iron_to_carbon
    POC_refr = d.particulate_refractory_fraction.carbon; PON_refr = d.particulate_refractory_fraction.nitrogen; POP_refr = d.particulate_refractory_fraction.phosphorus
    Qp(s,c,k)  = (s===:sp ? spP[c,k]/(spC[c,k]+εC)  : s===:diat ? diatP[c,k]/(diatC[c,k]+εC)  : diazP[c,k]/(diazC[c,k]+εC))
    Qfe(s,c,k) = (s===:sp ? spFe[c,k]/(spC[c,k]+εC) : s===:diat ? diatFe[c,k]/(diatC[c,k]+εC) : diazFe[c,k]/(diazC[c,k]+εC))
    results = CheckRow[]; C! = (a...; kw...) -> acompare!(results, a...; kw...)
    act = fin.(T); graze_active = act .& ((M["graze_sp"] .> 0) .| (M["graze_diat"] .> 0) .| (M["graze_diaz"] .> 0))
    C!("DON_prod = Q·f_toDON·DOC_prod", Q .* f_toDON .* M["DOC_prod"], M["DON_prod"]; mask = act)
    DOC_prod_pred = M["sp_loss_doc"] .+ M["diat_loss_doc"] .+ M["diaz_loss_doc"] .+ M["graze_sp_doc"] .+ M["graze_diat_doc"] .+ M["graze_diaz_doc"] .+ M["zoo_loss_doc"] .+ M["graze_zoo_doc"]
    C!("DOC_prod (routed-flux sum)", DOC_prod_pred, M["DOC_prod"]; mask = act)
    DOP_prod_pred = fill(NaN, ncols, nlev)
    for c in 1:ncols, k in 1:nlev
        act[c,k] || continue; Σdop = 0.0
        for (s, ag, al, aa, agp, alp, agz) in ((:sp,"graze_sp","sp_loss","sp_agg","graze_sp_poc","sp_loss_poc","graze_sp_zoo"),
                (:diat,"graze_diat","diat_loss","diat_agg","graze_diat_poc","diat_loss_poc","graze_diat_zoo"),
                (:diaz,"graze_diaz","diaz_loss","diaz_agg","graze_diaz_poc","diaz_loss_poc","graze_diaz_zoo"))
            q = Qp(s,c,k); pop = (M[agp][c,k] + M[alp][c,k] + M[aa][c,k]) * q
            R = (M[ag][c,k] + M[al][c,k] + M[aa][c,k]) * q - M[agz][c,k]*Qp_zoo - pop; Σdop += f_toDOP * max(R, 0.0)
        end
        DOP_prod_pred[c,k] = Qp_zoo*(M["zoo_loss_doc"][c,k] + M["graze_zoo_doc"][c,k]) + Σdop
    end
    C!("DOP_prod (§15.1 + §14.2)", DOP_prod_pred, M["DOP_prod"]; mask = graze_active)
    PFe_prod_pred = fill(NaN, ncols, nlev); PFe_prod_bio = fill(NaN, ncols, nlev)
    for c in 1:ncols, k in 1:nlev
        act[c,k] || continue; Σfe = 0.0
        for (s, aa, agp, alp) in ((:sp,"sp_agg","graze_sp_poc","sp_loss_poc"),(:diat,"diat_agg","graze_diat_poc","diat_loss_poc"),(:diaz,"diaz_agg","graze_diaz_poc","diaz_loss_poc"))
            Σfe += Qfe(s,c,k) * (M[aa][c,k] + M[agp][c,k] + M[alp][c,k])
        end
        PFe_prod_pred[c,k] = Σfe + Qfe_zoo * M["zoo_loss_poc"][c,k]
        PFe_prod_bio[c,k]  = M["P_iron_PROD"][c,k] - M["Fe_scavenge"][c,k] - Qfe_zoo * M["graze_zoo_poc"][c,k]
    end
    pfe_scale = maximum(x -> fin(x) ? abs(x) : 0.0, PFe_prod_bio); pfe_mask = graze_active .& (abs.(PFe_prod_bio) .> 1e-6 * pfe_scale)
    C!("PFe prod = P_iron_PROD − Fe_scav", PFe_prod_pred, PFe_prod_bio; mask = pfe_mask)
    light_cell = falses(ncols, nlev); dark_cell = falses(ncols, nlev)
    for c in 1:ncols
        φ = @view FRACR[:, c]; meanQSW = sum(φ .* @view QSW[:, c]); r = meanQSW > 0 ? (@view(QSW[:, c]) ./ meanQSW) : zero(φ)
        for k in 1:nlev
            act[c,k] || continue; lo = Inf; hi = -Inf
            for j in eachindex(φ); φ[j] > 0 || continue; PARj = PAR_avg[c,k]*r[j]; lo = min(lo, PARj); hi = max(hi, PARj); end
            lo > 1.1 && (light_cell[c,k] = true); hi < 0.9 && (dark_cell[c,k] = true)
        end
    end
    L = d.semilabile_remineralisation_rate_light; D = d.semilabile_remineralisation_rate_dark
    C!("DOC_remin light rate", DOC .* L.carbon,     M["DOC_remin"]; mask = light_cell)
    C!("DOC_remin dark rate",  DOC .* D.carbon,     M["DOC_remin"]; mask = dark_cell)
    C!("DON_remin light rate", DON .* L.nitrogen,   M["DON_remin"]; mask = light_cell)
    C!("DON_remin dark rate",  DON .* D.nitrogen,   M["DON_remin"]; mask = dark_cell)
    C!("DOP_remin light rate", DOP .* L.phosphorus, M["DOP_remin"]; mask = light_cell)
    C!("DOP_remin dark rate",  DOP .* D.phosphorus, M["DOP_remin"]; mask = dark_cell)
    R0 = d.refractory_remineralisation_rate; sub = falses(ncols, nlev); for c in 1:ncols, k in 2:nlev; sub[c,k] = act[c,k]; end
    C!("DOCr_remin (k≥2)", DOCr .* R0.carbon,     M["DOCr_remin"]; mask = sub)
    C!("DONr_remin (k≥2)", DONr .* R0.nitrogen,   M["DONr_remin"]; mask = sub)
    C!("DOPr_remin (k≥2)", DOPr .* R0.phosphorus, M["DOPr_remin"]; mask = sub)
    POC_tot = M["POC_REMIN_DIC"] .+ M["POC_REMIN_DOCr"]; POP_tot = M["POP_REMIN_PO4"] .+ M["POP_REMIN_DOPr"]; PON_tot = M["PON_REMIN_NH4"] .+ M["PON_REMIN_DONr"]
    C!("POC remin → DOCr fraction", fill(POC_refr, ncols, nlev), M["POC_REMIN_DOCr"] ./ POC_tot; mask = act .& (POC_tot .> 0))
    C!("POP remin → DOPr fraction", fill(POP_refr, ncols, nlev), M["POP_REMIN_DOPr"] ./ POP_tot; mask = act .& (POP_tot .> 0))
    C!("PON remin → DONr fraction", fill(PON_refr, ncols, nlev), M["PON_REMIN_DONr"] ./ PON_tot; mask = act .& (PON_tot .> 0))
    report_rows(results)
end

# =====================================================================================================
# Section 5 — carbonate & calcite (base config, sp the sole implicit calcifier)
# =====================================================================================================
function section_calcite()
    println("\n═══ 5. carbonate & calcite vs MARBL ═══")
    set!(tracers[:POC], flip(fill(0.0, ncols, nlev)))
    Msp_form = r2("sp_CaCO3_form"); Mtot_form = r2("CaCO3_form"); Mprod = r2("CaCO3_PROD"); MJspc = r2("J_spCaCO3")
    MH2CO3 = r2("H2CO3"); MCO3 = r2("CO3"); MpH = r2("pH_3D")
    DICm = r2("DIC"); ALKm = r2("ALK"); PO4m = r2("PO4"); SiO3m = r2("SiO3"); Tmm = r2("insitu_temp")
    pass = true
    println("── calcite formation / production (machine precision) ──")
    p = mreport("sp_CaCO3_form", mcompare((c,k,aux) -> calcite_formation(Val(:sp), c,1,k, grid, bgc, tracers, aux), Msp_form)...; tol = 1e-10); @test p; pass &= p
    p = mreport("CaCO3_form (Σ)", mcompare((c,k,aux) -> biological_calcium_carbonate_precipitation(c,1,k, grid, bgc.plankton, bgc, tracers, aux), Mtot_form)...; tol = 1e-10); @test p; pass &= p
    p = mreport("CaCO3_PROD", mcompare((c,k,aux) -> particulate_calcium_carbonate_production(c,1,k, grid, bgc.plankton, bgc, tracers, aux), Mprod)...; tol = 1e-10); @test p; pass &= p
    println("── living calcite tendency (machine precision) ──")
    p = mreport("J_spCaCO3", mcompare((c,k,aux) -> bgc(c,1,k, grid, Val(:spCaCO₃), clock, tracers, aux), MJspc)...; tol = 1e-8); @test p; pass &= p
    println("── carbonate speciation vs MARBL ──")
    marblP(c,lev) = Praw_ic[c,lev]
    println("   (1) CROSS-SOLVER — our full pH solve vs MARBL (tol = 1e-3 = MARBL's 3-s.f. limit)")
    for (klabel, ccx) in (("Perez&Fraga KF (default)", CC_PEREZ), ("Dickson&Riley KF (MARBL)", CC_DR))
        println("       KF = ", klabel)
        spk(c,lev) = carbon_dioxide_and_carbonate(ccx; DIC=DICm[c,lev], Alk=ALKm[c,lev], T=Tmm[c,lev], S=salinity[c,lev], phosphate=PO4m[c,lev], silicate=SiO3m[c,lev], P=marblP(c,lev))
        p = mreport("H2CO3 (CO₂aq)", mcompare((c,k,aux) -> first(spk(c, nlev - k + 1)), MH2CO3)...; tol = 1e-3); @test p; pass &= p
        p = mreport("CO3", mcompare((c,k,aux) -> last(spk(c, nlev - k + 1)), MCO3)...; tol = 1e-3); @test p; pass &= p
    end
    # (1b) the FULL MARBL-matched INDEPENDENT solve (CC_MARBL = DRTSafeSolver pH[6,9] + Dickson&Riley KF +
    # R=83.1451 + sulfate 96.062 + ρ=1026) reproduces MARBL's H2CO3/CO3 to MACHINE PRECISION — no shared H.
    # This is the whole point: an independent carbonate solve CAN match MARBL bit-for-bit (transcendental floor).
    println("   (1b) MARBL-MATCHED INDEPENDENT solve (CC_MARBL) vs MARBL — machine precision (tol 1e-9)")
    spkM(c,lev) = carbon_dioxide_and_carbonate(CC_MARBL; DIC=DICm[c,lev], Alk=ALKm[c,lev], T=Tmm[c,lev], S=salinity[c,lev], phosphate=PO4m[c,lev], silicate=SiO3m[c,lev], P=marblP(c,lev))
    p = mreport("H2CO3 (CC_MARBL)", mcompare((c,k,aux) -> first(spkM(c, nlev - k + 1)), MH2CO3)...; tol = 1e-9); @test p; pass &= p
    p = mreport("CO3 (CC_MARBL)",   mcompare((c,k,aux) -> last(spkM(c, nlev - k + 1)),  MCO3)...;  tol = 1e-9); @test p; pass &= p
    println("   (2) MACHINE PRECISION — our CO₂ & CO₃ ⟵ MARBL's equilibrium pH (tol 1e-10)")
    spH(c,lev) = carbon_dioxide_and_carbonate(CC_MARBL_R; DIC=DICm[c,lev], T=Tmm[c,lev], S=salinity[c,lev], pH=MpH[c,lev], phosphate=PO4m[c,lev], silicate=SiO3m[c,lev], P=marblP(c,lev))
    p = mreport("H2CO3 (CO₂aq)", mcompare((c,k,aux) -> first(spH(c, nlev - k + 1)), MH2CO3)...; tol = 1e-10); @test p; pass &= p
    p = mreport("CO3", mcompare((c,k,aux) -> last(spH(c, nlev - k + 1)), MCO3)...; tol = 1e-10); @test p; pass &= p
    println("   (3) MACHINE PRECISION — MARBL's CO₂ & CO₃ ⟶ our K₁K₂ equilibrium (tol 1e-10)")
    function co3_over_co2(c, lev)
        Tk = Tmm[c,lev] + 273.15; S = salinity[c,lev]; P = marblP(c,lev)
        K1 = CC_MARBL_R.carbonic_acid.K1(Tk, S; P); K2 = CC_MARBL_R.carbonic_acid.K2(Tk, S; P); H = 10.0 ^ (-MpH[c,lev]); K1 * K2 / H^2
    end
    p = mreport("CO3 ⟵ MARBL CO₂", mcompare((c,k,aux) -> (lev = nlev-k+1; MH2CO3[c,lev] * co3_over_co2(c, lev)), MCO3)...; tol = 1e-10); @test p; pass &= p
    p = mreport("CO₂ ⟵ MARBL CO3", mcompare((c,k,aux) -> (lev = nlev-k+1; MCO3[c,lev] / co3_over_co2(c, lev)), MH2CO3)...; tol = 1e-10); @test p; pass &= p
    pass
end

# =====================================================================================================
# Section 6 — O₂ / N redox (§18 nitrif/denitrif, §19 O₂). POC_REMIN_DIC supplied as the ballast input.
# =====================================================================================================
function section_oxygen()
    println("\n═══ 6. O₂ / N redox vs MARBL ═══")
    refract_C = bgc.detritus.particulate_refractory_fraction.carbon; pom_rate = bgc.detritus.particulate_organic_remineralisation_rate
    PDIC = r2("POC_REMIN_DIC")
    POC_inject = [fin(PDIC[c,l]) ? PDIC[c,l] / ((1 - refract_C) * pom_rate) : 0.0 for c in 1:ncols, l in 1:nlev]
    set!(tracers[:POC], flip(POC_inject))
    Mnit = r2("NITRIF"); Mden = r2("DENITRIF"); Mo2p = r2("O2_PRODUCTION"); Mo2c = r2("O2_CONSUMPTION")
    Mscf = r2("O2_CONSUMPTION_SCALEF"); MJo2 = r2("J_O2"); MJno3 = r2("J_NO3"); MJnh4 = r2("J_NH4")
    pass = true
    println("── rate diagnostics ──")
    p = mreport("O2_PRODUCTION", mcompare((c,k,aux) -> oxygen_production_total(c,1,k, grid, bgc, tracers, aux), Mo2p)...; tol = 1e-10); @test p; pass &= p
    p = mreport("NITRIF", mcompare((c,k,aux) -> water_column_nitrification(c,1,k, grid, bgc, tracers, aux), Mnit)...; tol = 1e-10); @test p; pass &= p
    p = mreport("DENITRIF (col1)", mcompare((c,k,aux) -> water_column_denitrification(c,1,k, grid, bgc, tracers, aux), Mden; colmask = c -> c == 1)...; tol = 1e-10); @test p; pass &= p
    Mo2c_scaled = [fin(Mo2c[c,lev]) && Mscf[c,lev] != 0 ? Mo2c[c,lev]/Mscf[c,lev] : FILL for c in 1:ncols, lev in 1:nlev]
    p = mreport("O2_CONS (/scalef)", mcompare((c,k,aux) -> oxygen_consumption(c,1,k, grid, bgc, tracers, aux), Mo2c_scaled)...; tol = 1e-10); @test p; pass &= p
    println("── assembled tendencies (real method) ──")
    p = mreport("J_O2 (scalef=1)", mcompare((c,k,aux) -> bgc(c,1,k, grid, Val(:O₂), clock, tracers, aux), MJo2;
                colmask = c -> all(l -> !active(c,l) || Mscf[c,l] == 1.0, 1:nlev))...; tol = 1e-8); @test p; pass &= p
    p = mreport("J_NO3", mcompare((c,k,aux) -> bgc(c,1,k, grid, Val(:NO₃), clock, tracers, aux), MJno3)...; tol = 1e-8); @test p; pass &= p
    p = mreport("J_NH4", mcompare((c,k,aux) -> bgc(c,1,k, grid, Val(:NH₄), clock, tracers, aux), MJnh4)...; tol = 1e-8); @test p; pass &= p
    pass
end

# =====================================================================================================
# Section 8 — iron cycle: 1-ligand speciation / scavenging / ligand budget (§17). ComplexedIron pure
# helpers fed MARBL's OWN state + sinking-flux diagnostics ⇒ machine precision, independent of our
# explicit-sinking flux approximation (which the operational model uses; the true ballast flux is Phase 11).
# =====================================================================================================
function section_iron()
    println("\n═══ 8. iron cycle (§17 speciation / scavenging / ligands) vs MARBL ═══")
    ci = ComplexedIron(); K = ci.ligand_stability_constant
    Fe = r2("Fe"); Lig = r2("Lig")
    Msc = Dict(n => r2(n) for n in
        ("Fefree","Fe_scavenge","Fe_scavenge_rate","Lig_scavenge","Lig_prod","Lig_deg","Lig_photochem","Lig_loss",
         "POC_sFLUX_IN","POC_hFLUX_IN","CaCO3_FLUX_IN","SiO2_FLUX_IN","dust_FLUX_IN",
         "POC_REMIN_DIC","POC_REMIN_DOCr","DOC_prod","photoFe_sp","photoFe_diat","photoFe_diaz"))
    keep = [active(c,lev) && lev != bottomlev(c) for c in 1:ncols, lev in 1:nlev]
    results = CheckRow[]; C! = (a...; kw...) -> acompare!(results, a...; kw...)

    # sinking mass (per area) reconstructed from MARBL's own POC/CaCO₃/bSi/dust flux_in diagnostics (§17.2 L2411).
    # ⚠️ `sinking_mass` takes **SI** fluxes (that is what the model feeds it), so MARBL's cgs diagnostics must be
    # converted here: mmol/m³·cm/s ⇒ ÷100, g/cm²/s ⇒ ×1e4. Feeding it MARBL's raw numbers "passes" while the
    # operational path is 100× out — a pure-helper test in the source's units cannot catch a units bug.
    Msink = fill(NaN, ncols, nlev)
    for c in 1:ncols, lev in 1:nlev
        keep[c,lev] || continue
        Fpoc = (Msc["POC_sFLUX_IN"][c,lev] + Msc["POC_hFLUX_IN"][c,lev]) / 100
        Msink[c,lev] = MP.sinking_mass(Fpoc, Msc["CaCO3_FLUX_IN"][c,lev] / 100,
                                       Msc["SiO2_FLUX_IN"][c,lev] / 100,
                                       Msc["dust_FLUX_IN"][c,lev] * 1e4, ci)
    end

    # §17.1 speciation + §17.2 scavenging (pure helpers, MARBL Fe/Lig + fluxes)
    Fefree_pred = fill(NaN, ncols, nlev); FeScav_pred = fill(NaN, ncols, nlev)
    FeScavRate_pred = fill(NaN, ncols, nlev); LigScav_pred = fill(NaN, ncols, nlev)
    for c in 1:ncols, lev in 1:nlev
        keep[c,lev] || continue
        Ff, FL = MP.iron_speciation(Fe[c,lev], Lig[c,lev], K)
        Fefree_pred[c,lev]     = Ff
        FeScavRate_pred[c,lev] = ci.iron_scavenging_rate * Msink[c,lev]
        FeScav_pred[c,lev]     = MP.iron_scavenging_flux(Ff, FL, Msink[c,lev], ci)
        LigScav_pred[c,lev]    = MP.ligand_scavenging_flux(FL, Msink[c,lev], ci)
    end
    C!("Fefree (§17.1 speciation)",  Fefree_pred,     Msc["Fefree"];           mask = keep)
    C!("Fe_scavenge_rate (r0·M)",    FeScavRate_pred, Msc["Fe_scavenge_rate"]; mask = keep)
    C!("Fe_scavenge (§17.2)",        FeScav_pred,     Msc["Fe_scavenge"];      mask = keep)
    C!("Lig_scavenge (§17.2)",       LigScav_pred,    Msc["Lig_scavenge"];     mask = keep)

    # §17.3 ligand production / bacterial degradation (pure, MARBL POC_remin + DOC_prod)
    POC_remin = Msc["POC_REMIN_DIC"] .+ Msc["POC_REMIN_DOCr"]
    LigProd_pred = ci.ligand_production_ratio .* (POC_remin .+ Msc["DOC_prod"])
    LigDeg_pred  = ci.ligand_degradation_rate .* POC_remin
    C!("Lig_prod (remin_to_Lig)",    LigProd_pred, Msc["Lig_prod"]; mask = keep)
    C!("Lig_deg (degrade_rate)",     LigDeg_pred,  Msc["Lig_deg"];  mask = keep)

    # §17.3 surface photo-oxidation: pure rate (interface-PAR reconstruction) × MARBL Lig
    LigPhoto_pred = fill(NaN, ncols, nlev)
    for c in 1:ncols, k in 1:nlev
        lev = nlev - k + 1; keep[c,lev] || continue
        LigPhoto_pred[c,lev] = Lig[c,lev] * MP.ligand_photochemistry_rate(c, 1, k, grid, auxs[c])
    end
    C!("Lig_photochem (surface UV)", LigPhoto_pred, Msc["Lig_photochem"]; mask = keep)

    # §17.3 total ligand loss = Lig_scav + 0.20·Σ photoFe + Lig_photochem + Lig_deg (validates the 0.20 factor + assembly)
    LigLoss_pred = LigScav_pred .+ 0.20 .* (Msc["photoFe_sp"] .+ Msc["photoFe_diat"] .+ Msc["photoFe_diaz"]) .+ LigPhoto_pred .+ LigDeg_pred
    C!("Lig_loss (§17.3 assembly)",  LigLoss_pred, Msc["Lig_loss"]; mask = keep)

    # §20.1 assembled ligand tendency J_Lig = Lig_prod − Lig_loss (from the machine-precision components)
    C!("J_Lig (= Lig_prod − Lig_loss)", LigProd_pred .- LigLoss_pred, r2("J_Lig"); mask = keep)
    rows_pass = report_rows(results)

    # §20.1 assembled dissolved-iron tendency J_Fe (real method) vs MARBL. J_Fe = P_iron_remin − Σ photoFe −
    # Fe_scavenge + auto/zoo DIC/DOC routing. The ballast P_iron_remin is provisional (Phase 11), so inject
    # PFe with PFe·pom_rate = MARBL's P_iron_REMIN (as §6 does for POC), and add MARBL's Fe_scavenge back to
    # the target (our operational scavenging uses the inferred sinking flux). This validates the ENTIRE
    # dissolved-Fe routing (photoFe + the §14.2 auto/zoo Fe→DIC/DOC/POC split + remin) end-to-end on the
    # plain-Fe model (scavenging off ⇒ ∂Fe = the biological J_Fe).
    println("── assembled dissolved-iron tendency (real method) ──")
    pom_rate = bgc.detritus.particulate_organic_remineralisation_rate
    PIR = r2("P_iron_REMIN")
    set!(tracers[:PFe], flip([fin(PIR[c,l]) ? PIR[c,l] / pom_rate : 0.0 for c in 1:ncols, l in 1:nlev]))
    JFe_biology = r2("J_Fe") .+ r2("Fe_scavenge")
    pj = mreport("J_Fe (biology; +Fe_scav, PFe inj)",
                 mcompare((c,k,aux) -> bgc(c,1,k, grid, Val(:Fe), clock, tracers, aux), JFe_biology)...; tol = 1e-8)
    @test pj
    return rows_pass & pj
end

# =====================================================================================================
# Section 7 — +cocco config (Eppley temp / VCO2 / picpoc), freshly-built cocco baseline + own model
# =====================================================================================================
function section_cocco()
    println("\n═══ 7. +cocco config (Eppley / VCO2 / picpoc) vs MARBL ═══")
    HISTc = joinpath(MARBL_ROOT, "MARBL", "tests", "input_files", "baselines", "call_compute_subroutines.cocco.history.nc")
    isfile(HISTc) || (println("  [skip] +cocco baseline not present ($HISTc) — build it to enable this section."); return true)
    dsc = NCDataset(HISTc)
    r2c(n) = permutedims(Array{Float64}(coalesce.(dsc[n][:, :], FILL)))
    zt_c = Array(dsc["zt"][:]) / 100; zw_c = Array(dsc["zw"][:]) / 100
    nc = size(r2c("insitu_temp"), 1); nl = length(zt_c)
    gridc = RectilinearGrid(size = (nc, 1, nl), x = (0, nc), y = (0, 1), z = vcat(-reverse(zw_c), 0.0), topology = (Periodic, Periodic, Bounded))
    flipc(A) = (B = fill(0.0, nc, 1, nl); for c in 1:nc, k in 1:nl; v = A[c, nl-k+1]; B[c,1,k] = fin(v) ? v : 0.0; end; B)
    PARf = CenterField(gridc); Sf = CenterField(gridc)
    bgc0c = NutrientsPlanktonDetritus(gridc;
            nutrients = Nutrients(NitrateAmmonia(; nitrification_rate = 0.0), PO₄, Fe, Si),
            plankton  = MARBLPlankton_cocco(gridc), detritus = MARBLDetritus(gridc; sinking_speed = 10.0),
            inorganic_carbon = ExplicitCalciumCarbonate(gridc; calcium_carbonate_dissolution_rate = 0.03/day, carbon_chemistry = CONST_RHO_CC),
            oxygen = MARBLOxygen(), light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PARf))
    modelc = NonhydrostaticModel(gridc; tracers = PHYSICS_TRACERS, biogeochemistry = bgc0c, buoyancy = nothing, auxiliary_fields = (S = Sf,))
    bgcc = modelc.biogeochemistry.underlying_biogeochemistry
    marbl_of_c = Dict(:NO₃=>"NO3", :NH₄=>"NH4", :PO₄=>"PO4", :Fe=>"Fe", :Si=>"SiO3", :DIC=>"DIC", :Alk=>"ALK", :O₂=>"O2", :T=>"insitu_temp")
    zerofill_c = Set((:POP, :PFe, :bSi, :POC, :CaCO₃)); trc = modelc.tracers
    for name in (required_biogeochemical_tracers(bgcc)..., PHYSICS_TRACERS...); name in zerofill_c && continue; set!(trc[name], flipc(r2c(get(marbl_of_c, name, ascii_name(name))))); end
    set!(PARf, flipc(r2c("PAR_avg")))
    set!(bgcc.plankton.carbon_dioxide, flipc(r2c("H2CO3")))    # prescribe MARBL's CO₂aq (VCO2 + picpoc input)
    clockc = modelc.clock
    subc(c) = (φ = @view FRACR[:, c]; mq = sum(φ .* @view QSW[:, c]); mq > 0 ? ntuple(j -> (φ[j], (@view(QSW[:, c]) ./ mq)[j]), nbin) : ntuple(j -> (j == 1 ? 1.0 : 0.0, 0.0), nbin))
    PARif = build_PARiface(r2c("spChl") .+ r2c("diatChl") .+ r2c("diazChl") .+ r2c("coccoChl"), r2c("PAR_avg"), gridc)
    auxsc = ntuple(c -> (PAR = SubcolumnPAR(PARf, map(first, subc(c)), map(last, subc(c))),
                         PAR_interface = PARif), nc)
    Tmc = r2c("insitu_temp"); activec(c,lev) = fin(Tmc[c,lev]); botc(c) = maximum(lev for lev in 1:nl if activec(c,lev))
    function ccompare(pred, truth; floor_frac = 1e-9)
        keep(c,lev) = activec(c,lev) && lev != botc(c) && fin(truth[c,lev])
        scale = 0.0; for c in 1:nc, lev in 1:nl; keep(c,lev) && (scale = max(scale, abs(truth[c,lev]))); end
        thr = floor_frac * scale; mr = 0.0; n = 0; worst = (0,0,0.0,0.0)
        for c in 1:nc, k in 1:nl
            lev = nl - k + 1; keep(c,lev) || continue; abs(truth[c,lev]) ≤ thr && continue
            p = pred(c, k, auxsc[c]); isfinite(p) || continue; n += 1
            r = abs(p - truth[c,lev]) / (abs(truth[c,lev]) + 1e-14); r > mr && (mr = r; worst = (c,lev,p,truth[c,lev]))
        end
        (mr, n, worst)
    end
    pass = true
    println("── cocco growth + calcification (machine precision) ──")
    p = mreport("photoC_cocco", ccompare((c,k,aux) -> photosynthesis(Val(:cocco), c,1,k, gridc, bgcc.plankton, trc, aux), r2c("photoC_cocco"))...; tol = 1e-10); @test p; pass &= p
    p = mreport("cocco_CaCO3_form", ccompare((c,k,aux) -> calcite_formation(Val(:cocco), c,1,k, gridc, bgcc, trc, aux), r2c("cocco_CaCO3_form"))...; tol = 1e-10); @test p; pass &= p
    println("── assembled cocco element tendencies (machine precision) ──")
    for (nm, tr) in (("J_coccoC", r2c("J_coccoC")), ("J_coccoChl", r2c("J_coccoChl")), ("J_coccoP", r2c("J_coccoP")), ("J_coccoFe", r2c("J_coccoFe")), ("J_coccoCaCO3", r2c("J_coccoCaCO3")))
        sym = tracer_symbol(nm[3:end])
        p = mreport(nm, ccompare((c,k,aux) -> bgcc(c,1,k, gridc, Val(sym), clockc, trc, aux), tr)...; tol = 1e-8); @test p; pass &= p
    end

    # ---- sediment at the +cocco sea floor: the sediment is config-independent (a pure function of the floor
    # ---- fluxes), so this checks it plugs into the cocco tracer set and closes every bottom-cell J. cocco (not
    # ---- sp) is the calcifier, so CaCO₃ burial is fed by the cocco-sourced pool. Same reconstruction as base.
    println("── +cocco assembled bottom-cell tendencies: bgc(Val) + sediment/Δz  vs  MARBL J_<tracer> ──")
    r1c(n) = Array{Float64}(coalesce.(dsc[n][:], FILL))
    Δzm_c = [zw_bot - zw_top for (zw_top, zw_bot) in zip(vcat(0.0, zw_c)[1:nl], zw_c)]
    sedc = MARBLSediment()
    Sc = reconstruct_sediment_inputs(sedc, r2c, r1c, botc, Δzm_c, zw_c .* 100, nc, bgcc)
    pass &= assembled_bottom_cell_J(sedc, Sc, r2c, n -> haskey(dsc, n), gridc, bgcc, trc, clockc, botc, nc, nl,
                                    auxsc, COCCO_J_TRACERS)

    # ---- +cocco implicit ballast: the sweep must reproduce every cocco-config flux/remin/tendency too ----
    pass &= ballast_config("+cocco", dsc, r2c, flipc, gridc, nc, nl, bgcc.plankton, PARf, Sf,
                           marbl_of_c, auxsc, activec, botc; prescribe_co2 = true)
    close(dsc)
    pass
end

# =====================================================================================================
# Section 9 — sediments: burial / sedimentary denitrification / "other" remin (§16.5, Phase 9).
# Every check is at the sea floor (lev = bottomlev(c) = MARBL's kmt), the one cell every other section
# masks out. Our sediment is a pure function of the floor flux, so feeding it MARBL's OWN `pocToFloor` /
# `calcToFloor` / bottom O₂,NO₃,Ω isolates the §16.5 formulas from our explicit-sinking flux (Phase 11).
# =====================================================================================================

# absolute error scaled by the largest |truth| across the five columns (so an exact-zero truth — e.g.
# calcToSed in an undersaturated column — is still a meaningful check). Machine precision ⇔ ≤ 1e-13.
function sreport(name, ours, truth; tol = 1e-12)
    scale = maximum(abs, truth)
    mr = maximum(abs(o - t) for (o, t) in zip(ours, truth)) / (scale + 1e-300)
    pass = mr ≤ tol
    tag = mr ≤ 1e-13 ? "MACHINE-PRECISION ✅" : pass ? "PASS ✅" : "FAIL ❌"
    @printf("  %-22s max-rel=%.2e over %d columns  %s\n", name, mr, length(truth), tag)
    pass || println("      ours =", ours, "\n      marbl=", truth)
    @test pass
    pass
end

# Ragueneau (2000) bSi burial length scale + O₂ scale factor (marbl_interior_tendency_mod L2960-2990),
# and MARBL's dust→Fe molar conversion — hoisted so both the base and +cocco sediment checks share them.
const SCALELEN_Z = [100.0e2, 250.0e2, 500.0e2, 1000.0e2]
const SCALELEN_V = [1.0, 3.6, 4.7, 4.8]
function scalelen(zwk)
    zwk <  SCALELEN_Z[1]   && return SCALELEN_V[1]
    zwk >= SCALELEN_Z[end] && return SCALELEN_V[end]
    for n in 2:4
        zwk < SCALELEN_Z[n] && return SCALELEN_V[n-1] +
            (SCALELEN_V[n]-SCALELEN_V[n-1]) * (zwk-SCALELEN_Z[n-1]) / (SCALELEN_Z[n]-SCALELEN_Z[n-1])
    end
end
o2_scalefac(o2) = o2 < 45.0 ? 1 + (2.6 - 1) * min(1.0, (45.0 - o2) / (45.0 - 5.0)) : 1.0
const DUST_TO_FE = 0.035 / 55.845 * 1.0e9

# Reconstruct MARBL's five particulate floor fluxes at the sea floor from its own flux diagnostics: POC/CaCO₃
# are diagnosed (`pocToFloor`/`calcToFloor`); PFe/POP/bSi are reconstructed (CESM2.1 has POP%gamma = 0,
# parm_SiO2_gamma = 0 ⇒ purely soft), plus the bottom-cell O₂/NO₃/Ω and a `sedflux(c, name)` closure that runs
# the discrete-form sediment on those inputs. Config-independent — shared by the base and +cocco sections.
function reconstruct_sediment_inputs(sed, r2x, r1x, botx, Δzm_x, zw_cm_x, ncx, bgcx)
    cm2 = 100.0
    dzcm(c) = Δzm_x[botx(c)] * 100
    pocFloor = r1x("pocToFloor"); pocSed = r1x("pocToSed"); calcFloor = r1x("calcToFloor")
    popSed = r1x("popToSed"); bsiSed = r1x("bsiToSed"); calcSed = r1x("calcToSed")
    O2a = r2x("O2"); NO3a = r2x("NO3"); CO3a = r2x("CO3"); CSa = r2x("co3_sat_calc")
    O2b  = [O2a[c, botx(c)]  for c in 1:ncx]
    NO3b = [NO3a[c, botx(c)] for c in 1:ncx]
    Ωb   = [CO3a[c, botx(c)] / CSa[c, botx(c)] for c in 1:ncx]   # MARBL's test: CO3 > thres·CO3_sat_calcite

    Pfi = r2x("P_iron_FLUX_IN"); Ppr = r2x("P_iron_PROD"); Prm = r2x("P_iron_REMIN")
    dRem = r2x("dust_REMIN"); fesed = r2x("FESEDFLUX")
    pfeFloor = [(l = botx(c); Pfi[c,l] + dzcm(c) * (Ppr[c,l] - (Prm[c,l] - dRem[c,l]*DUST_TO_FE - fesed[c,l]/dzcm(c))))
                for c in 1:ncx]

    refract_C = bgcx.detritus.particulate_refractory_fraction.carbon
    SFI = r2x("SiO2_FLUX_IN"); SPR = r2x("SiO2_PROD")
    PPI = r2x("POP_FLUX_IN"); PPR = r2x("POP_PROD"); CFI = r2x("POC_FLUX_IN"); PDIC = r2x("POC_REMIN_DIC")
    bsiFloor = [(l = botx(c); ℓ = scalelen(zw_cm_x[l]) * o2_scalefac(O2b[c]) * 650e2; d = exp(-dzcm(c)/ℓ);
                 SFI[c,l]*d + SPR[c,l]*(1-d)*ℓ) for c in 1:ncx]
    popFloor = [(l = botx(c);
                 poc_remin_int = PDIC[c,l]/(1 - refract_C) - (pocFloor[c] - pocSed[c])/dzcm(c);
                 PPI[c,l] + dzcm(c) * (PPR[c,l] - poc_remin_int * PPI[c,l]/CFI[c,l])) for c in 1:ncx]

    # O₂_consumption_scale = 1 here: the sediment Val{:O₂} returns UNSCALED consumption, and `assembled` below
    # applies MARBL's scalef itself (j = p − (p−j)·Mscf), so it must not be double-scaled inside the sediment.
    tracked(c) = (; O₂ = fill(O2b[c]), NO₃ = fill(NO3b[c]), Ω = fill(Ωb[c]), O₂_consumption_scale = fill(1.0),
                    POC = fill(pocFloor[c]/cm2), POP = fill(popFloor[c]/cm2), bSi = fill(bsiFloor[c]/cm2),
                    CaCO₃ = fill(calcFloor[c]/cm2), PFe = fill(pfeFloor[c]/cm2))
    sedflux(c, name) = sed(1, 1, nothing, Val(name), nothing, nothing, tracked(c))

    (; dzcm, pocFloor, pocSed, popFloor, popSed, bsiFloor, bsiSed, calcFloor, calcSed, pfeFloor,
       O2b, NO3b, Ωb, refract_C, sedflux)
end

# The end-to-end Phase-9 check: our pointwise `bgc(Val(X))` PLUS the sediment coupling (÷Δz, exactly as
# `update_tendencies!` adds it) vs MARBL's `J_<X>` at the sea floor. Each sinking pool is injected so our
# first-order interior remin equals MARBL's INTERIOR remin (its *_REMIN diagnostics already contain the
# sediment return at kmt, so we back that out), leaving the return flux to come from the sediment.
function assembled_bottom_cell_J(sed, S, r2x, hasx, gridx, bgcx, tracersx, clockx, botx, ncx, nlx, auxsx, jlist)
    cm2 = 100.0
    dzcm = S.dzcm
    refract_P = bgcx.detritus.particulate_refractory_fraction.phosphorus
    pom = bgcx.detritus.particulate_organic_remineralisation_rate
    si  = bgcx.detritus.particulate_silica_remineralisation_rate
    ca  = bgcx.inorganic_carbon.calcium_carbonate_dissolution_rate[1]   # one rate per realisation
    PDIC = r2x("POC_REMIN_DIC"); POP_PO4 = r2x("POP_REMIN_PO4")
    SiO2_REM = r2x("SiO2_REMIN"); CaCO3_REM = r2x("CaCO3_REMIN"); PFe_REM = r2x("P_iron_REMIN")

    flipx(A) = (B = fill(0.0, ncx, 1, nlx);
                for c in 1:ncx, k in 1:nlx; v = A[c, nlx-k+1]; B[c,1,k] = fin(v) ? v : 0.0; end; B)
    function inject!(name, total_of, rate, returned)
        A = fill(0.0, ncx, nlx)
        for c in 1:ncx, l in 1:nlx
            t = total_of(c, l); fin(t) || continue
            A[c, l] = ((l == botx(c)) ? t - returned[c]/dzcm(c) : t) / rate
        end
        set!(tracersx[name], flipx(A))
    end
    inject!(:POC,   (c,l) -> PDIC[c,l]/(1 - S.refract_C), pom, S.pocFloor .- S.pocSed)
    inject!(:POP,   (c,l) -> POP_PO4[c,l]/(1 - refract_P), pom, S.popFloor .- S.popSed)
    inject!(:bSi,   (c,l) -> SiO2_REM[c,l],  si, S.bsiFloor .- S.bsiSed)
    inject!(:CaCO₃, (c,l) -> CaCO3_REM[c,l], ca, S.calcFloor .- S.calcSed)
    inject!(:PFe,   (c,l) -> PFe_REM[c,l],   pom, zeros(ncx))

    # ExplicitCalciumCarbonate's DIC/Alk methods read an Ω auxiliary field (abiotic precipitation, rate 0 here so its
    # value is immaterial — but the field must exist); the water column isn't run so supply a zero field
    Ωfield = CenterField(gridx)
    auxsΩ = ntuple(c -> merge(auxsx[c], (Ω = Ωfield,)), ncx)

    Mscf = hasx("O2_CONSUMPTION_SCALEF") ? r2x("O2_CONSUMPTION_SCALEF") : nothing   # absent ⇒ scalef = 1
    Fescav = r2x("Fe_scavenge")
    coupled = MP.coupled_tracers(sed)
    function assembled(c, sym)
        l = botx(c); k = nlx - l + 1
        j = bgcx(c, 1, k, gridx, Val(sym), clockx, tracersx, auxsΩ[c])
        sym in coupled && (j += S.sedflux(c, sym) * cm2 / dzcm(c))
        # MARBL scales the whole O₂ consumption (interior + sediment) by O2_CONSUMPTION_SCALEF, production is not
        sym === :O₂ && Mscf !== nothing && (p = oxygen_production_total(c, 1, k, gridx, bgcx, tracersx, auxsΩ[c]);
                                            j = p - (p - j) * Mscf[c, l])
        sym === :Fe && (j -= Fescav[c, l])                 # plain Fe: scavenging is not a bgc sink
        j
    end

    # MARBL's 10-day NO₃-depletion clamp (compute_denitrif L3405-3409) scales down denitrif+sed_denitrif TOGETHER
    # where NO₃ < 10·spd·(denitrif+sed_denitrif). Our split sediment/oxygen model clamps only the water-column
    # denitrif, not sed_denitrif (a documented Phase-9 deviation), so at a near-NO₃-depleted floor cell the three
    # denitrification-coupled tendencies (NO₃/Alk/O₂) differ. Detect those cells from MARBL's own diagnostics and
    # exclude them for those tracers only (everything else, including sed_denitrif itself, stays machine precise).
    spd = 86400.0
    DEN = r2x("DENITRIF")
    clamp_bites(c) = (l = botx(c);
        sd_vol = sedimentary_denitrification(sed, S.pocFloor[c]/cm2, S.O2b[c], S.NO3b[c]) / (dzcm(c)/100);
        S.NO3b[c] < 10 * spd * (DEN[c,l] + sd_vol))
    clamped = filter(clamp_bites, 1:ncx)
    isempty(clamped) || println("     NO₃-clamp active at column(s) ", join(clamped, ","),
                                " — NO₃/Alk/O₂ excluded there (unimplemented sed_denitrif clamp, §16.5 deviation)")

    pass = true
    for (sym, jn) in jlist
        hasx(jn) || continue
        J = r2x(jn)
        cols = sym in (:NO₃, :Alk, :O₂) ? [c for c in 1:ncx if !clamp_bites(c)] : collect(1:ncx)
        isempty(cols) && continue
        pass &= sreport(jn, [assembled(c, sym) for c in cols], [J[c, botx(c)] for c in cols]; tol = 1e-9)
    end
    pass
end

function section_sediments()
    println("\n═══ 9. sediments (§16.5 burial / sed. denitrification / other remin) vs MARBL ═══")
    sed = MARBLSediment()
    r1(n) = Array{Float64}(coalesce.(ds[n][:], FILL))
    cm2 = 100.0                                   # mmol m⁻² s⁻¹ → MARBL's nmol cm⁻² s⁻¹
    kmt(c) = bottomlev(c)

    S = reconstruct_sediment_inputs(sed, r2, r1, kmt, Δz_m, zw .* 100, ncols, bgc)
    Fpoc(c) = S.pocFloor[c] / cm2

    # every burial term is validated the SAME way: feed our burial fraction MARBL's own floor flux (diagnosed
    # for POC/CaCO₃, reconstructed for PFe/POP/bSi) and reproduce MARBL's `*ToSed` diagnostic to machine
    # precision. `buried_X` is the tendency, an area flux (mmol m⁻² s⁻¹ → MARBL's nmol cm⁻² s⁻¹ via ×100).
    pass = true
    println("── burial fluxes to the sediment (fed MARBL's own floor fluxes) ──")
    pass &= sreport("pocToSed", [S.sedflux(c, :buried_C)     * cm2 for c in 1:ncols], S.pocSed)
    pass &= sreport("ponToSed", [S.sedflux(c, :buried_N)     * cm2 for c in 1:ncols], r1("ponToSed"))
    pass &= sreport("popToSed", [S.sedflux(c, :buried_P)     * cm2 for c in 1:ncols], S.popSed)
    pass &= sreport("bsiToSed", [S.sedflux(c, :buried_Si)    * cm2 for c in 1:ncols], S.bsiSed)
    pass &= sreport("pfeToSed", [S.sedflux(c, :buried_Fe)    * cm2 for c in 1:ncols], r1("pfeToSed"))
    # Ω = 0.894 (col 1, buried) vs 0.776 (col 2, dissolved) straddles caco3_bury_thres_omega_calc = 0.89:
    # this fails outright if the formulation-doc value of 1.0 is used instead of the source value
    println("     Ω at the floor: ", join((@sprintf("%.3f", Ω) for Ω in S.Ωb), ", "), "  (threshold ", sed.calcite_burial_saturation_threshold, ")")
    pass &= sreport("calcToSed", [S.sedflux(c, :buried_CaCO₃) * cm2 for c in 1:ncols], S.calcSed)

    println("── sedimentary denitrification & anaerobic remineralisation ──")
    pass &= sreport("SedDenitrif",
                    [sedimentary_denitrification(sed, Fpoc(c), S.O2b[c], S.NO3b[c]) * cm2 for c in 1:ncols],
                    r1("SedDenitrif"))
    pass &= sreport("OtherRemin",
                    [other_remineralisation(sed, Fpoc(c), S.O2b[c], S.NO3b[c]) * cm2 for c in 1:ncols],
                    r1("OtherRemin"))

    # ω = 0 at every floor cell (bottom O₂ ≥ 117 mmol/m³) ⇒ no water-column denitrification response;
    # sed_denitrif enters NO₃ once, UNSCALED by ω (MARBL J_NO3 = nitrif − denitrif − sed_denitrif − uptake)
    println("── NO₃: sedimentary denitrification enters once, unscaled ──")
    Mden = r2("DENITRIF")
    pass &= sreport("DENITRIF(kmt)", [Mden[c, kmt(c)] for c in 1:ncols], zeros(ncols); tol = 0.0)
    pass &= sreport("denitrif response",
                    [denitrification_response(sed, Fpoc(c), S.O2b[c], S.NO3b[c]) for c in 1:ncols],
                    zeros(ncols); tol = 0.0)
    pass &= sreport("our ∂NO₃|sed", [S.sedflux(c, :NO₃) * cm2 for c in 1:ncols], -r1("SedDenitrif"))

    # ---- ASSEMBLED bottom-cell tendencies vs MARBL J_<tracer> at kmt, for EVERY prognostic tracer of the base
    # ---- config. The 8 sediment-coupled tracers get the return flux; the DOM + all 15 PFT pools receive none,
    # ---- so they double as a leak guard. `J_Lig`/`_ALT_CO2` are ComplexedIron / the 2nd realization (not this
    # ---- config); POC/POP/bSi/CaCO₃/PFe have no MARBL `J_` (implicit ballast flux, our sinking tracers).
    println("── assembled bottom-cell tendencies: bgc(Val) + sediment/Δz  vs  MARBL J_<tracer> ──")
    pass &= assembled_bottom_cell_J(sed, S, r2, n -> haskey(ds, n), grid, bgc, tracers, clock, kmt, ncols, nlev,
                                    auxs, BASE_J_TRACERS)

    pass
end

# every prognostic tracer of the base config, paired with its MARBL J_ diagnostic (sp calcifies → spCaCO3)
const BASE_J_TRACERS =
    ((:DIC,"J_DIC"), (:Alk,"J_ALK"), (:NH₄,"J_NH4"), (:PO₄,"J_PO4"), (:Si,"J_SiO3"),
     (:NO₃,"J_NO3"), (:O₂,"J_O2"), (:Fe,"J_Fe"),
     (:DOCr,"J_DOCr"), (:DONr,"J_DONr"), (:DOPr,"J_DOPr"), (:DOC,"J_DOC"), (:DON,"J_DON"), (:DOP,"J_DOP"),
     (:spC,"J_spC"), (:spChl,"J_spChl"), (:spFe,"J_spFe"), (:spCaCO₃,"J_spCaCO3"), (:spP,"J_spP"),
     (:diatC,"J_diatC"), (:diatChl,"J_diatChl"), (:diatFe,"J_diatFe"), (:diatP,"J_diatP"), (:diatSi,"J_diatSi"),
     (:diazC,"J_diazC"), (:diazChl,"J_diazChl"), (:diazFe,"J_diazFe"), (:diazP,"J_diazP"),
     (:zooC,"J_zooC"))

# +cocco config: sp is a NON-calcifier (no spCaCO3), the calcite pool is carried by coccolithophores
const COCCO_J_TRACERS =
    ((:DIC,"J_DIC"), (:Alk,"J_ALK"), (:NH₄,"J_NH4"), (:PO₄,"J_PO4"), (:Si,"J_SiO3"),
     (:NO₃,"J_NO3"), (:O₂,"J_O2"), (:Fe,"J_Fe"),
     (:DOCr,"J_DOCr"), (:DONr,"J_DONr"), (:DOPr,"J_DOPr"), (:DOC,"J_DOC"), (:DON,"J_DON"), (:DOP,"J_DOP"),
     (:spC,"J_spC"), (:spChl,"J_spChl"), (:spFe,"J_spFe"), (:spP,"J_spP"),
     (:diatC,"J_diatC"), (:diatChl,"J_diatChl"), (:diatFe,"J_diatFe"), (:diatP,"J_diatP"), (:diatSi,"J_diatSi"),
     (:diazC,"J_diazC"), (:diazChl,"J_diazChl"), (:diazFe,"J_diazFe"), (:diazP,"J_diazP"),
     (:coccoC,"J_coccoC"), (:coccoChl,"J_coccoChl"), (:coccoFe,"J_coccoFe"), (:coccoP,"J_coccoP"),
     (:coccoCaCO₃,"J_coccoCaCO3"), (:zooC,"J_zooC"))

# =====================================================================================================
# Section 10 — the GENERAL N_A/N_Z configs (Phase 10) vs MARBL runs of the SAME configs.
#
# Both baselines are built by general_config/make_baselines.jl, which generates MARBL's settings file from
# the very same Julia structs the model here is built from (so the two sides cannot drift by a typo), and
# runs the shipped `marbl.exe`. `general5` turns on everything Phase 10 adds — 5 autotrophs (two calcifiers,
# two silicifiers), 2 zooplankton, MULTI-MEMBER prey classes, ZOO-ON-ZOO grazing, and prey eaten by both
# predators; `general2` is the degenerate end (2 autotrophs, 1 predator, no N-fixer, no carbon limitation).
#
# What is compared: every per-PFT rate, the WHOLE grazing matrix (per prey: total grazing, the split to each
# predator separately, and the POC/DOC routes; per predator: assimilation, loss and the food-weighted
# f_zoo_detr), and the assembled tendency of every PLANKTON tracer. The plankton tracers are the ones Phase
# 10 rewrote and they close the grazing matrix; the aggregate DIC/POC/DOC/nutrient tendencies additionally
# need MARBL's implicit ballast remineralisation (Phase 11), which the base-config sections above supply by
# injecting POC_REMIN_DIC — that gap is a property of Phase 11, not of the N_A/N_Z generality, so it is not
# re-litigated here.
# =====================================================================================================

function section_general(name, plankton_of, asnames, zsnames; prescribe_co2)
    println("\n═══ 10. general config `$name` ($(length(asnames)) autotrophs, $(length(zsnames)) zooplankton) vs MARBL ═══")
    HIST = joinpath(DRAFT_DIR, "general_config", "baselines", "$name.history.nc")
    isfile(HIST) || (println("  [skip] baseline missing ($HIST) — run general_config/make_baselines.jl"); return true)
    dsg = NCDataset(HIST)
    r2g(n) = permutedims(Array{Float64}(coalesce.(dsg[n][:, :], FILL)))
    zt_g = Array(dsg["zt"][:]) / 100; zw_g = Array(dsg["zw"][:]) / 100
    ng = size(r2g("insitu_temp"), 1); nlg = length(zt_g)
    gridg = RectilinearGrid(size = (ng, 1, nlg), x = (0, ng), y = (0, 1),
                            z = vcat(-reverse(zw_g), 0.0), topology = (Periodic, Periodic, Bounded))
    flipg(A) = (B = fill(0.0, ng, 1, nlg); for c in 1:ng, k in 1:nlg; v = A[c, nlg-k+1]; B[c,1,k] = fin(v) ? v : 0.0; end; B)

    PARg = CenterField(gridg); Sg = CenterField(gridg)
    bgc0g = NutrientsPlanktonDetritus(gridg;
            nutrients = Nutrients(NitrateAmmonia(; nitrification_rate = 0.0), PO₄, Fe, Si),
            plankton  = plankton_of(gridg), detritus = MARBLDetritus(gridg; sinking_speed = 10.0),
            inorganic_carbon = ExplicitCalciumCarbonate(gridg; calcium_carbonate_dissolution_rate = 0.03/day, carbon_chemistry = CONST_RHO_CC),
            oxygen = MARBLOxygen(), light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PARg))
    modelg = NonhydrostaticModel(gridg; tracers = PHYSICS_TRACERS, biogeochemistry = bgc0g, buoyancy = nothing, auxiliary_fields = (S = Sg,))
    bgcg = modelg.biogeochemistry.underlying_biogeochemistry
    pg = bgcg.plankton; trg = modelg.tracers; clockg = modelg.clock

    marbl_of_g = Dict(:NO₃=>"NO3", :NH₄=>"NH4", :PO₄=>"PO4", :Fe=>"Fe", :Si=>"SiO3",
                      :DIC=>"DIC", :Alk=>"ALK", :O₂=>"O2", :T=>"insitu_temp")

    # MARBL's `autotroph_zero_consistency_enforce` (marbl_interior_tendency_mod.F90:1048, called from
    # setup_local_tracers:793 — i.e. BEFORE PAR and every rate): if any one of an autotroph's Chl/C/P/Fe
    # (/Si, if it silicifies) is exactly zero in a cell, ALL of that autotroph's tracers are zeroed there.
    # The shipped IC has such cells at depth (a nonzero C beside an exactly-zero Chl), so MARBL's rates see
    # a dead autotroph where the raw state has a live one. MARBL's tendency is a function of the ENFORCED
    # state, so that is the state to feed in. (The base sections instead drop those cells via
    # `domcleanmask`; enforcing is the stronger check — it keeps every cell in the comparison.)
    enforced = Dict{Symbol, Matrix{Float64}}()
    for s in asnames
        flags = MP.traits(getproperty(pg.autotrophs, s))
        mask = falses(ng, nlg)
        for nm in (MP.carbon_name(s), MP.chlorophyll_name(s), MP.phosphorus_name(s), MP.iron_name(s))
            mask .|= r2g(ascii_name(nm)) .== 0
        end
        flags.silicifier && (mask .|= r2g(ascii_name(MP.silicon_name(s))) .== 0)
        for nm in MP.autotroph_tracer_names(s, flags)
            A = r2g(ascii_name(nm)); A[mask] .= 0.0; enforced[nm] = A
        end
    end
    stateg(nm) = get(enforced, nm, nothing) === nothing ? r2g(get(marbl_of_g, nm, ascii_name(nm))) : enforced[nm]

    zerofill_g = Set((:POP, :PFe, :bSi, :POC, :CaCO₃))
    for nm in (required_biogeochemical_tracers(bgcg)..., PHYSICS_TRACERS...)
        nm in zerofill_g && continue
        set!(trg[nm], flipg(stateg(nm)))
    end
    set!(PARg, flipg(r2g("PAR_avg")))
    prescribe_co2 && set!(pg.carbon_dioxide, flipg(r2g("H2CO3")))   # MARBL's CO₂aq (cocco VCO2 + picpoc)

    subg(c) = (φ = @view FRACR[:, c]; mq = sum(φ .* @view QSW[:, c]);
               mq > 0 ? ntuple(j -> (φ[j], (@view(QSW[:, c]) ./ mq)[j]), nbin) :
                        ntuple(j -> (j == 1 ? 1.0 : 0.0, 0.0), nbin))
    totChl = sum(stateg(MP.chlorophyll_name(s)) for s in asnames)
    PARifg = build_PARiface(totChl, r2g("PAR_avg"), gridg)
    auxg = ntuple(c -> (PAR = SubcolumnPAR(PARg, map(first, subg(c)), map(last, subg(c))),
                        PAR_interface = PARifg), ng)

    Tmg = r2g("insitu_temp"); activeg(c,lev) = fin(Tmg[c,lev]); botg(c) = maximum(lev for lev in 1:nlg if activeg(c,lev))
    function gcompare(pred, truth; floor_frac = 1e-9)
        keep(c,lev) = activeg(c,lev) && lev != botg(c) && fin(truth[c,lev])
        scale = 0.0; for c in 1:ng, lev in 1:nlg; keep(c,lev) && (scale = max(scale, abs(truth[c,lev]))); end
        thr = floor_frac * scale; mr = 0.0; n = 0; worst = (0,0,0.0,0.0)
        for c in 1:ng, k in 1:nlg
            lev = nlg - k + 1; keep(c,lev) || continue; abs(truth[c,lev]) ≤ thr && continue
            p = pred(c, k, auxg[c]); isfinite(p) || continue; n += 1
            r = abs(p - truth[c,lev]) / (abs(truth[c,lev]) + 1e-14); r > mr && (mr = r; worst = (c,lev,p,truth[c,lev]))
        end
        (mr, n, worst)
    end
    check(nm, f, truth; tol = 1e-9) = (p = mreport(nm, gcompare(f, truth)...; tol); @test p; p)

    pass = true
    println("── per-autotroph growth, calcification & loss ──")
    for s in asnames
        V = Val(s)
        pass &= check("photoC_$s", (c,k,a) -> photosynthesis(V, c,1,k, gridg, pg, trg, a), r2g("photoC_$s"))
        pass &= check("$(s)_loss", (c,k,a) -> MP.autotroph_linear_loss(V, c,1,k, gridg, pg, trg), r2g("$(s)_loss"))
        pass &= check("$(s)_agg", (c,k,a) -> MP.autotroph_aggregation(V, c,1,k, gridg, pg, trg), r2g("$(s)_agg"))
        pass &= check("$(s)_loss_poc", (c,k,a) -> MP.loss_to_poc(V, c,1,k, gridg, pg, trg), r2g("$(s)_loss_poc"))
        haskey(dsg, "$(s)_CaCO3_form") &&
            (pass &= check("$(s)_CaCO3_form", (c,k,a) -> calcite_formation(V, c,1,k, gridg, bgcg, trg, a),
                           r2g("$(s)_CaCO3_form")))
    end

    # the grazing matrix, prey by prey — autotrophs AND zooplankton. `graze_<prey>_<pred>` is MARBL's
    # per-predator split, so a prey eaten by two predators (and a prey class whose members share one flux by
    # biomass) is checked term by term, not just in the sum.
    println("── grazing matrix: per-prey totals, per-PREDATOR split, and POC/DOC routing ──")
    for m in (asnames..., zsnames...)
        V = Val(m)
        pass &= check("graze_$m", (c,k,a) -> MP.grazing_loss(V, c,1,k, gridg, pg, trg), r2g("graze_$m"))
        pass &= check("graze_$(m)_zootot", (c,k,a) -> MP.graze_to_zooplankton(V, c,1,k, gridg, pg, trg), r2g("graze_$(m)_zootot"))
        pass &= check("graze_$(m)_poc", (c,k,a) -> MP.graze_to_poc(V, c,1,k, gridg, pg, trg), r2g("graze_$(m)_poc"))
        pass &= check("graze_$(m)_doc", (c,k,a) -> MP.graze_to_doc(V, c,1,k, gridg, pg, trg), r2g("graze_$(m)_doc"))
        for z in zsnames
            pass &= check("graze_$(m)_$(z)", (c,k,a) -> MP.graze_to_predator(V, Val(z), c,1,k, gridg, pg, trg),
                          r2g("graze_$(m)_$(z)"))
        end
    end

    println("── per-predator assimilation, loss & food-weighted detrital routing ──")
    for z in zsnames
        V = Val(z)
        pass &= check("x_graze_$z", (c,k,a) -> MP.zooplankton_grazing_gain(V, c,1,k, gridg, pg, trg), r2g("x_graze_$z"))
        pass &= check("$(z)_loss", (c,k,a) -> MP.zooplankton_loss(V, c,1,k, gridg, pg, trg), r2g("$(z)_loss"))
        pass &= check("$(z)_loss_poc", (c,k,a) -> MP.zoo_loss_to_poc(V, c,1,k, gridg, pg, trg), r2g("$(z)_loss_poc"))
        pass &= check("$(z)_loss_doc", (c,k,a) -> MP.zoo_loss_to_doc(V, c,1,k, gridg, pg, trg), r2g("$(z)_loss_doc"))
    end

    println("── assembled plankton tendencies: bgc(Val(:x)) vs MARBL J_<tracer> ──")
    plankton_tracers = filter(nm -> nm !== :T, required_biogeochemical_tracers(pg))
    for nm in plankton_tracers
        pass &= check("J_" * ascii_name(nm), (c,k,a) -> bgcg(c,1,k, gridg, Val(nm), clockg, trg, a), r2g("J_" * ascii_name(nm)); tol = 1e-8)
    end

    # ---- implicit ballast for this general config: sweep flux/remin + every assembled tendency vs MARBL ----
    pass &= ballast_config(name, dsg, r2g, flipg, gridg, ng, nlg, pg, PARg, Sg,
                           marbl_of_g, auxg, activeg, botg; prescribe_co2 = prescribe_co2)

    close(dsg)
    return pass
end

# =====================================================================================================
# Reusable implicit-ballast comparison — ONE `MARBLBallastDetritus` column driven by a config's MARBL
# state, checked against that config's Fortran flux/remin/tendency diagnostics AND for element closure.
# Called for base CESM, +cocco, general2 (2P1Z) and general5 (5P2Z) so every config exercises the sweep.
#
# Units (MARBL is cgs):  FLUX_IN "mmol/m^3 cm/s" ⇒ ÷100 for mmol/m²/s;  dust "g/cm^2/s" ⇒ ×1e4 for g/m²/s;
#                        REMIN "mmol/m^3/s" is already SI.
# Index/sign: MARBL's FLUX_IN(lev) is the (positive-down) flux INTO cell `lev` through its TOP face. Our
# cell is k = nlev-lev+1, whose top face is k+1, and our F ≤ 0 downward ⇒ magnitude = -F[c, 1, k+1].
#
# `open_bottom = false`: the floor flux is remineralised in the bottom cell (nothing buried), so the column
# CLOSES — which is what makes the element-conservation check below exact. This changes ONLY the bottom
# cell, and every comparison excludes it (`lev != botX`), so the interior match to MARBL is untouched.
function ballast_config(tag, dsX, r2X, flipX, gridX, ncX, nlX, plankton, PARfieldX, SfieldX,
                        marbl_of_X, auxX, activeX, botX; prescribe_co2 = false)
    println("\n── [$tag] implicit ballast: fluxes + remin + 𝓜 + every assembled tendency vs MARBL J_<tracer> ──")

    # DUST_FLUX drives the QA-dust deficit 𝒬 (⇒ the hard-POC ballast); FESEDFLUX is the sedimentary iron
    # source added to P_iron%remin below the clip. Both are per-cm² in MARBL (×1e4 / ÷100 to SI).
    DF    = Array{Float64}(coalesce.(dsX["DUST_FLUX"][:], FILL))                 # g/cm²/s (surface)
    dustF = Field{Center, Center, Nothing}(gridX)
    for c in 1:ncX; @inbounds dustF[c, 1, 1] = fin(DF[c]) ? DF[c] * 1e4 : 0.0; end
    FSF   = r2X("FESEDFLUX")
    fesed = CenterField(gridX)
    set!(fesed, flipX([fin(FSF[c, l]) ? FSF[c, l] / 100 : 0.0 for c in 1:ncX, l in 1:nlX]))

    bd   = MARBLBallastDetritus(gridX;
                                ballast = MARBLBallast(; surface_dust_flux = dustF,
                                                         sedimentary_iron_flux = fesed, open_bottom = false),
                                flux_diagnostics = true)
    bgcB = NutrientsPlanktonDetritus(gridX;
            nutrients = Nutrients(NitrateAmmonia(; nitrification_rate = 0.0), PO₄, ComplexedIron(), Si),
            plankton  = plankton, detritus = bd,
            inorganic_carbon = ImplicitExplicitCalcite(gridX; carbon_chemistry = CONST_RHO_CC),
            oxygen = MARBLOxygen(),
            light_attenuation = PrescribedPhotosyntheticallyActiveRadiation(PARfieldX))
    modelB = NonhydrostaticModel(gridX; tracers = PHYSICS_TRACERS, biogeochemistry = bgcB, buoyancy = nothing,
                                 auxiliary_fields = (S = SfieldX,))
    bgcBu = modelB.biogeochemistry.underlying_biogeochemistry
    trB   = modelB.tracers; clkB = modelB.clock
    for name in (required_biogeochemical_tracers(bgcBu)..., PHYSICS_TRACERS...)
        set!(trB[name], flipX(r2X(get(marbl_of_X, name, ascii_name(name)))))
    end
    prescribe_co2 && set!(bgcBu.plankton.carbon_dioxide, flipX(r2X("H2CO3")))
    # fire ONLY the detritus sweep — NOT the full fan-out. The plankton fan-out step recomputes the aqueous
    # CO₂ from carbon chemistry, which would clobber MARBL's prescribed H2CO3 (and, with S = 0 in the +cocco
    # column, with a wrong value); the calcite/CO₂-coupled tendencies then diverge. The sweep uses the
    # prescribed CO₂ for Π_CaCO₃, and nothing here needs the ImplicitExplicitCalcite Ω (sediment-free interior).
    update_biogeochemical_state!(modelB, bgcBu.detritus, bgcBu)      # ← fires compute_marbl_ballast!

    F = bd.fluxes; RM = bd.remineralisation; Rc = bgcBu.inorganic_carbon.remineralisation

    # config-local compare engine (same shape as ccompare/gcompare): every active interior cell, machine
    # precision, with an optional per-column mask for the (unimplemented) O₂ consumption scale factor.
    function xcmp(pred, truth; colmask = c -> true, floor_frac = 1e-9)
        keep(c, lev) = colmask(c) && activeX(c, lev) && lev != botX(c) && fin(truth[c, lev])
        scale = 0.0; for c in 1:ncX, lev in 1:nlX; keep(c, lev) && (scale = max(scale, abs(truth[c, lev]))); end
        thr = floor_frac * scale; mr = 0.0; n = 0; worst = (0, 0, 0.0, 0.0)
        for c in 1:ncX, k in 1:nlX
            lev = nlX - k + 1; keep(c, lev) || continue; abs(truth[c, lev]) ≤ thr && continue
            p = pred(c, k, auxX[c]); isfinite(p) || continue; n += 1
            r = abs(p - truth[c, lev]) / (abs(truth[c, lev]) + 1e-14); r > mr && (mr = r; worst = (c, lev, p, truth[c, lev]))
        end
        (mr, n, worst)
    end
    chk(nm, f, truth; tol = 1e-9, colmask = c -> true) = (p = mreport(nm, xcmp(f, truth; colmask)...; tol); @test p; p)
    fx(fld, sc)  = (c, k, a) -> -fld[c, 1, k + 1] * sc     # flux magnitude at the cell's top face
    cellf(fld)   = (c, k, a) -> fld[c, 1, k]

    pass = true
    println("   fluxes (mmol/m³·cm/s; dust g/cm²/s)")
    pass &= chk("POC_sFLUX_IN",   fx(F.POC_soft, 100.0), r2X("POC_sFLUX_IN"))
    pass &= chk("POC_hFLUX_IN",   fx(F.POC_hard, 100.0), r2X("POC_hFLUX_IN"))
    pass &= chk("POC_FLUX_IN",    (c, k, a) -> -(F.POC_soft[c, 1, k + 1] + F.POC_hard[c, 1, k + 1]) * 100.0, r2X("POC_FLUX_IN"))
    pass &= chk("CaCO3_FLUX_IN",  fx(F.CaCO₃, 100.0),    r2X("CaCO3_FLUX_IN"))
    pass &= chk("SiO2_FLUX_IN",   fx(F.bSi,   100.0),    r2X("SiO2_FLUX_IN"))
    pass &= chk("dust_FLUX_IN",   fx(F.dust,  1e-4),     r2X("dust_FLUX_IN"))
    pass &= chk("P_iron_FLUX_IN", fx(F.PFe,   100.0),    r2X("P_iron_FLUX_IN"))
    pass &= chk("POP_FLUX_IN",    fx(F.POP,   100.0),    r2X("POP_FLUX_IN"))
    println("   sinking mass 𝓜 (= Fe_scavenge_rate / parm_Fe_scavenge_rate0)")
    pass &= chk("sinking_mass",   cellf(bd.sinking_mass), r2X("Fe_scavenge_rate") ./ ComplexedIron().iron_scavenging_rate)
    println("   remin (mmol/m³/s)")
    pass &= chk("POC_remin",   cellf(RM.POC), r2X("POC_REMIN_DIC") .+ r2X("POC_REMIN_DOCr"))
    pass &= chk("POP_remin",   cellf(RM.POP), r2X("POP_REMIN_PO4") .+ r2X("POP_REMIN_DOPr"))
    pass &= chk("SiO2_REMIN",  cellf(RM.bSi), r2X("SiO2_REMIN"))
    pass &= chk("P_iron_REMIN",cellf(RM.PFe), r2X("P_iron_REMIN"))
    pass &= chk("CaCO3_REMIN", cellf(Rc),     r2X("CaCO3_REMIN"))

    # the strongest gate: the ACTUAL assembled tendency the timestepper integrates, for EVERY prognostic
    # tracer, vs MARBL's J_<tracer> — proves the ballast remin is routed correctly into every reservoir.
    println("   assembled tendency bgc(Val:x) vs MARBL J_<tracer>")
    o2col = c -> true
    if haskey(dsX, "O2_CONSUMPTION_SCALEF")                # base only; o2_consumption_scalef unimplemented
        Mscf = r2X("O2_CONSUMPTION_SCALEF")
        o2col = c -> all(l -> !activeX(c, l) || Mscf[c, l] == 1.0, 1:nlX)
    end
    for nm in (required_biogeochemical_tracers(bgcBu)..., PHYSICS_TRACERS...)
        nm === :T && continue
        truth = r2X("J_" * get(marbl_of_X, nm, ascii_name(nm)))
        cm = nm === :O₂ ? o2col : (c -> true)
        pass &= chk("J_$nm", (c, k, a) -> bgcBu(c, 1, k, gridX, Val(nm), clkB, trB, a), truth; tol = 1e-8, colmask = cm)
    end

    # NB conservation is NOT checked here. The element-summed column integral on THIS grid does not vanish —
    # iron has genuine external sources (dust dissolution + fesedflux), and where MARBL's sea floor sits above
    # the grid bottom the sinking flux escapes into below-floor cells — and matching MARBL's own `Jint_<E>tot`
    # budget would require replicating its sediment burial (which has no counterpart in this sediment-free
    # ballast column). The routing IS proven to conserve, transitively: every assembled `J_<tracer>` above
    # matches MARBL (which conserves) to machine precision. Direct C/N/P/Fe/Si closure is checked on a properly
    # closed box, per plankton config, in `test_phase11.jl` (the established box pattern).

    # `P_CaCO3_ALT_CO2%prod = P_CaCO3%prod` (MARBL L2524) ⇒ calcite flux/remin are provably identical across
    # replicates, which is why ONE remin field serves them all; assert MARBL's own two agree.
    alt = true
    for c in 1:ncX, lev in 1:nlX
        (activeX(c, lev) && lev != botX(c) && fin(r2X("CaCO3_FLUX_IN")[c, lev])) || continue
        alt &= r2X("CaCO3_ALT_CO2_FLUX_IN")[c, lev] == r2X("CaCO3_FLUX_IN")[c, lev]
        alt &= r2X("CaCO3_ALT_CO2_REMIN")[c, lev]   == r2X("CaCO3_REMIN")[c, lev]
    end
    @printf("   %-12s %s\n", "ALT_CO2≡CaCO3", alt ? "IDENTICAL ✅" : "DIFFER ❌"); @test alt
    return pass & alt
end

# =====================================================================================================
# Section 11 — implicit ballast sinking (Armstrong 2000): the flux profiles and remin ARE the model
# =====================================================================================================
# §8 reconstructs the sinking mass from MARBL's OWN `POC_sFLUX_IN`/`CaCO3_FLUX_IN`/… diagnostics, so it
# validates the scavenging *formula* while taking the fluxes on trust. This section is what earns them:
# a `MARBLBallastDetritus` column with NO particulate tracers, whose sweep must reproduce every one of
# MARBL's flux and remin profiles from the tracer state alone.
#
# Units (MARBL is cgs):  FLUX_IN "mmol/m^3 cm/s" ⇒ ÷100 for mmol/m²/s;  dust "g/cm^2/s" ⇒ ×1e4 for g/m²/s;
#                        REMIN "mmol/m^3/s" is already SI.
# Index/sign: MARBL's FLUX_IN(lev) is the (positive-down) flux INTO cell `lev` through its TOP face. Our
# cell is k = nlev-lev+1, whose top face is k+1, and our F ≤ 0 downward ⇒ magnitude = -F[c, 1, k+1].
function section_ballast()
    println("\n═══ 11. implicit ballast sinking (Armstrong 2000) — flux + remin profiles vs MARBL ═══")
    return ballast_config("base CESM", ds, r2, flip, grid, ncols, nlev, MARBLPlankton(),
                          PARfield, Sfield, marbl_of, auxs, active, bottomlev)
end

# =====================================================================================================
# run all sections
# =====================================================================================================
println("MARBL numerical comparison — OceanBioME vs standalone MARBL Fortran driver")
println("(base CESM2.1 60-level column; ice-radiation sub-columns + interface PAR prescribed)")

ALLPASS = Ref(true)
@testset "MARBL comparison vs Fortran driver (all phases)" begin
    ALLPASS[] &= section_rates()
    ALLPASS[] &= section_realmethod()
    ALLPASS[] &= section_jtend()
    ALLPASS[] &= section_dompom()
    ALLPASS[] &= section_calcite()
    ALLPASS[] &= section_oxygen()
    ALLPASS[] &= section_iron()
    ALLPASS[] &= section_sediments()
    ALLPASS[] &= section_ballast()
    close(ds)
    ALLPASS[] &= section_cocco()
    ALLPASS[] &= section_general("general2", (g) -> general2_plankton(),
                                 (:sp, :diat), (:zoo, ); prescribe_co2 = false)
    ALLPASS[] &= section_general("general5", general5_plankton,
                                 (:sp, :diat, :diaz, :cocco, :diat2), (:zoo, :zoo2); prescribe_co2 = true)
end

println()
println(ALLPASS[] ? "ALL MARBL comparison sections PASS ✅" : "SOME MARBL comparison checks FAILED ❌")
ALLPASS[] || error("MARBL comparison failed")

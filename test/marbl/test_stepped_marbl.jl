##### STEPPED FULL-MODEL GATE — the strongest MARBL check: build each MARBL column on its OWN kmt-deep grid,
##### step the WHOLE model through the real Oceananigans timestepper (update_state! → compute_tendencies! →
##### update_tendencies! sediment coupling ⇒ model.timestepper.Gⁿ), and assert Gⁿ reproduces MARBL's
##### interior_tendencies `J_<tracer>` to MACHINE PRECISION at every interior cell — with the sediment coupling
##### injected at the true sea floor. advection=nothing + no closure ⇒ Gⁿ is pure biogeochemical source/sink
##### (+ sediment flux/Δz at the floor), directly comparable to MARBL's interior tendencies. Run for ALL FOUR
##### configs: base CESM (sp/diat/diaz), +cocco (4P1Z), general2 (2P1Z), general5 (5P2Z).
#####
##### Per-column kmt-deep grids make sediment.bottom_indices = 1 = the true sea floor, so the implicit ballast
##### sweep runs surface→floor through only physical cells (no FILL-poisoned below-floor cells) and the sediment
##### deposits its remineralised return at the correct cell. Exact MARBL light is supplied through the model's
##### own light path (cell-mean PAR + interface PAR + ice sub-columns). Carbonate uses the independent CC_MARBL
##### solve (drtsafe pH[6,9] + Dickson&Riley KF + R=83.1451); O₂ uses MARBL's O2_CONSUMPTION_SCALEF.
#####
##### Run: julia --check-bounds=yes --project=marbl-info3/draft/cmpenv marbl-info3/draft/test_stepped_marbl.jl

using NCDatasets, Printf, Test, Oceananigans, OceanBioME
using Oceananigans.Fields: CenterField, Field, set!, compute!, AbstractField
using Oceananigans.Grids: Center, Face
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_transition
import Oceananigans.Biogeochemistry: biogeochemical_auxiliary_fields, update_biogeochemical_state!
using Oceananigans.TimeSteppers: update_state!



include(joinpath(@__DIR__, "marbl_names.jl"))
include(joinpath(@__DIR__, "marbl_baselines.jl"))
const MP = @__MODULE__
require_baselines("test_stepped_marbl") || exit(0)

include(joinpath(@__DIR__, "general_configs.jl"))    # general2_plankton / general5_plankton
using OceanBioME: NutrientsPlanktonDetritus, Nutrients, NitrateAmmonia, PO₄, Si,
                  PrescribedPhotosyntheticallyActiveRadiation, CarbonChemistry
const CCM = OceanBioME.Models.CarbonChemistryModel

# ---------------------------------------------------------------------------------------------------
# baselines + shared column physics (all 4 configs are the SAME 5 columns; salinity + ice bins are shared)
# ---------------------------------------------------------------------------------------------------
const MARBLDIR = joinpath(MARBL_ROOT, "MARBL", "tests", "input_files")
const HIST_base  = joinpath(MARBLDIR, "baselines", "call_compute_subroutines.history.nc")
const HIST_cocco = joinpath(MARBLDIR, "baselines", "call_compute_subroutines.cocco.history.nc")
const HIST_gen(name) = joinpath(DRAFT_DIR, "general_config", "baselines", "$name.history.nc")
const ICF  = joinpath(MARBLDIR, "initial_conditions", "call_compute_subroutines.20190718.nc")
const FILL = 9.96920996838687e36
fin(x) = (x != FILL) & isfinite(x)
_median(v) = (s = sort(collect(v)); n = length(s); iseven(n) ? (s[n÷2] + s[n÷2+1]) / 2 : s[n÷2+1])

# shared physics from the base IC file: salinity + ice-radiation sub-columns (nbin, col)
let d = NCDataset(ICF)
    Sraw = permutedims(Array{Float64}(coalesce.(d["salinity"][:, :], FILL)))
    global salinity = (10 < _median(filter(fin, Sraw)) < 45) ? Sraw : Sraw .* 1000.0
    global Praw = permutedims(Array{Float64}(coalesce.(d["pressure"][:, :], FILL)))   # MARBL's prescribed pressure (bar)
    global FRACR = Array{Float64}(d["FRACR_BIN"][:, :]); global QSW = Array{Float64}(d["QSW_BIN"][:, :]); close(d)
end
const nbin = size(FRACR, 1)
subcols_of(c) = (φ = @view FRACR[:, c]; meanQSW = sum(φ .* @view QSW[:, c]);
                 meanQSW > 0 ? ntuple(j -> (φ[j], (@view(QSW[:, c]) ./ meanQSW)[j]), nbin) :
                               ntuple(j -> (j == 1 ? 1.0 : 0.0, 0.0), nbin))

# the sub columns MARBL's ice radiation bins describe: one φⱼ and one rⱼ per bin, each of
# which may be a scalar (uniform over the domain — every MARBL water column is its own grid/model here) or a
# per-(i, j) field. φⱼ/rⱼ are the ice-radiation area fraction and shortwave ratio, PARⱼ(k) = PAR_avg(k)·rⱼ.
# the ice radiation bins, as sub column weights carried on the light itself
subcolumn_par(field, c) = (s = subcols_of(c);
                           SubcolumnPAR(field, ntuple(j -> s[j][1], nbin), ntuple(j -> s[j][2], nbin)))
KPARdz(chl, dz_m) = (w = max(chl, 0.02); (w < 0.13224 ? 0.000919 * w^0.3536 : 0.001131 * w^0.4562) * (dz_m * 100))

marbl_of = Dict(:NO₃=>"NO3", :NH₄=>"NH4", :PO₄=>"PO4", :Fe=>"Fe", :Si=>"SiO3", :DIC=>"DIC", :Alk=>"ALK",
                :O₂=>"O2", :T=>"insitu_temp")

# ---- MARBL-matched carbon chemistry (see project_carbonate_solver_comparison) ----
struct DicksonRileyKF{PC}; pressure_correction :: PC; end
@inline (c::DicksonRileyKF)(T, S, Is, KS; P = nothing) =
    c.pressure_correction(T, P) *
    exp(1590.2 / T - 12.641 + 1.525 * sqrt(Is) + log(1 - 0.001005 * S) + log(1 + (0.14 / 96.062 * S / 1.80655) / KS))
const CC_MARBL = CarbonChemistry(; density_function = (args...) -> 1026.0,
    solver = OceanBioME.DRTSafeSolver(; max_iters = 100, bracket_grows = 3, atol = 1e-10, H_lower = 1e-9, H_upper = 1e-6),
    fluoride = DicksonRileyKF(CCM.PressureCorrection(; a₀=-9.78, a₁=-0.0090, a₂=-0.000942, b₀=-0.00391, b₁=0.000054)),
    carbonic_acid = (
        K1 = CCM.K1(; pressure_correction = CCM.PressureCorrection(; a₀=-25.50, a₁=0.1271, a₂=0.0, b₀=-0.00308, b₁=0.0000877, R=83.1451)),
        K2 = CCM.K2(; pressure_correction = CCM.PressureCorrection(; a₀=-15.82, a₁=-0.0219, a₂=0.0, b₀=0.00113, b₁=-0.0001475, R=83.1451))))

# ---- MARBL's exact light, prescribed through the model's own light path
# (aux.PAR / PAR_interface / PAR_subcolumns) ----
struct StepLight{NT}; aux :: NT; end
biogeochemical_auxiliary_fields(l::StepLight) = l.aux
update_biogeochemical_state!(model, l::StepLight) = (for f in values(l.aux); f isa AbstractField && compute!(f); end; nothing)

# MARBLPlankton*/general*_plankton @eval their per-PFT tendency methods on construction (manifest_*!); do it ONCE
# at top level so per-column models (built inside stepped_config) call them at a valid world age, not silently zero.
const DUMMY = RectilinearGrid(size = (1, 1, 1), extent = (1, 1, 1))
MARBLPlankton(); MARBLPlankton_cocco(DUMMY); general2_plankton(); general5_plankton(DUMMY)

# ---------------------------------------------------------------------------------------------------
# step one config: build each column on its own kmt-deep grid, update_state!, compare Gⁿ to MARBL J
# ---------------------------------------------------------------------------------------------------
function stepped_config(tag, HISTX, plankton_of, asnames)
    isfile(HISTX) || (@warn "baseline missing — skipping $tag ($HISTX)"; return)
    dsX = NCDataset(HISTX)
    r2X(n) = permutedims(Array{Float64}(coalesce.(dsX[n][:, :], FILL)))
    zwX = Array(dsX["zw"][:]) / 100
    TmX = r2X("insitu_temp"); ncX = size(TmX, 1); nlX = length(zwX)
    activeX(c, l) = fin(TmX[c, l]); botX(c) = maximum(l for l in 1:nlX if activeX(c, l))
    Δz_mX = [b - t for (t, b) in zip(vcat(0.0, zwX)[1:nlX], zwX)]
    O2SCF = haskey(dsX, "O2_CONSUMPTION_SCALEF") ? r2X("O2_CONSUMPTION_SCALEF") : nothing
    DFs   = Array{Float64}(coalesce.(dsX["DUST_FLUX"][:], FILL))
    FSFs  = r2X("FESEDFLUX")
    PARavg = r2X("PAR_avg")
    totChl = sum(r2X(ascii_name(MP.chlorophyll_name(s))) for s in asnames)

    # MARBL autotroph_zero_consistency_enforce: if any of an autotroph's C/Chl/P/Fe(/Si) is 0 in a cell, ALL of
    # that autotroph's tracers are zeroed there (marbl_interior_tendency_mod.F90:1048). The shipped IC has such
    # cells; MARBL's J is a function of the ENFORCED state, so feed that in.
    plankton0 = plankton_of(RectilinearGrid(size = (1, 1, nlX), x = (0, 1), y = (0, 1), z = vcat(-reverse(zwX), 0.0), topology = (Periodic, Periodic, Bounded)))
    enforced = Dict{Symbol, Matrix{Float64}}()
    for s in asnames
        flags = MP.traits(getproperty(plankton0.autotrophs, s))
        mask = falses(ncX, nlX)
        for nm in (MP.carbon_name(s), MP.chlorophyll_name(s), MP.phosphorus_name(s), MP.iron_name(s))
            mask .|= r2X(ascii_name(nm)) .== 0
        end
        flags.silicifier && (mask .|= r2X(ascii_name(MP.silicon_name(s))) .== 0)
        for nm in MP.autotroph_tracer_names(s, flags)
            A = r2X(ascii_name(nm)); A[mask] .= 0.0; enforced[nm] = A
        end
    end
    stateX(nm) = get(enforced, nm, nothing) === nothing ? r2X(get(marbl_of, nm, ascii_name(nm))) : enforced[nm]

    function pariface_col(c, kmt, gc)
        Pif = Field{Center, Center, Face}(gc); Iface = fill(0.0, kmt + 1)
        if fin(PARavg[c,1]) && PARavg[c,1] > 0
            Iface[1] = PARavg[c,1] * KPARdz(totChl[c,1], Δz_mX[1]) / (1 - exp(-KPARdz(totChl[c,1], Δz_mX[1])))
            for l in 1:kmt; Iface[l+1] = Iface[l] * exp(-KPARdz(totChl[c,l], Δz_mX[l])); end
        end
        for kf in 1:(kmt + 1); @inbounds Pif[1,1,kf] = Iface[kmt - kf + 2]; end
        Pif
    end

    istat   = Dict{Symbol, Vector{Float64}}()
    min_inj = Dict(nm => Inf for nm in (:O₂,:NO₃,:PO₄,:Si,:Fe,:DIC,:Alk))
    maxint_assembly = 0.0; max_pocfloor = 0.0
    println("▶ $tag: $ncX columns, building/stepping..."); flush(stdout)

    for c in 1:ncX
        kmt = botX(c)
        gc  = RectilinearGrid(size = (1, 1, kmt), x = (0, 1), y = (0, 1),
                              z = vcat(-reverse(zwX[1:kmt]), 0.0), topology = (Periodic, Periodic, Bounded))
        cflip(v) = reshape([v[kmt - k + 1] for k in 1:kmt], 1, 1, kmt)
        lev_of(k) = kmt - k + 1

        PARc  = CenterField(gc); set!(PARc, cflip([PARavg[c,l]   for l in 1:kmt]))
        Sc    = CenterField(gc); set!(Sc,   cflip([salinity[c,l] for l in 1:kmt]))
        dustc = Field{Center, Center, Nothing}(gc); @inbounds dustc[1,1,1] = fin(DFs[c]) ? DFs[c] * 1e4 : 0.0
        fesc  = CenterField(gc); set!(fesc, cflip([fin(FSFs[c,l]) ? FSFs[c,l] / 100 : 0.0 for l in 1:kmt]))
        o2scf = CenterField(gc); set!(o2scf, cflip([(O2SCF !== nothing && fin(O2SCF[c,l])) ? O2SCF[c,l] : 1.0 for l in 1:kmt]))
        Pc    = CenterField(gc); set!(Pc, cflip([Praw[c,l] for l in 1:kmt]))   # MARBL's prescribed pressure for the carbonate
        lgt   = StepLight((PAR = subcolumn_par(PARc, c), PAR_interface = pariface_col(c, kmt, gc)))

        bdc  = MARBLBallastDetritus(gc; ballast = MARBLBallast(; surface_dust_flux = dustc, sedimentary_iron_flux = fesc, open_bottom = true))
        bgcc = NutrientsPlanktonDetritus(gc;
                  nutrients = Nutrients(NitrateAmmonia(; nitrification_rate = 0.0), PO₄, ComplexedIron(), Si),
                  plankton  = plankton_of(gc), detritus = bdc,
                  inorganic_carbon = ImplicitExplicitCalcite(gc; carbon_chemistry = CC_MARBL),
                  oxygen = MARBLOxygen(; oxygen_consumption_scale_factor = o2scf),
                  sediment = MARBLSedimentModel(gc; implicit_sinking = true),
                  light_attenuation = lgt)
        mc = NonhydrostaticModel(gc; tracers = PHYSICS_TRACERS, biogeochemistry = bgcc, advection = nothing, buoyancy = nothing, auxiliary_fields = (S = Sc, pressure = Pc))
        for nm in (required_biogeochemical_tracers(mc.biogeochemistry.underlying_biogeochemistry)...,
                   PHYSICS_TRACERS...)
            set!(mc.tracers[nm], cflip([stateX(nm)[c,l] for l in 1:kmt]))
        end

        update_state!(mc)
        Gc = mc.timestepper.Gⁿ
        fc = Oceananigans.fields(mc)
        sm = mc.biogeochemistry.sediment

        # compare EVERY active cell k = 1..kmt. The floor cell (k=1 ⇔ MARBL lev kmt) carries the sediment
        # coupling in Gⁿ, and MARBL's J_<tracer> there includes the benthic flux too — so it must match as well.
        for nm in keys(Gc)
            nm === :T && continue
            jn = "J_" * get(marbl_of, nm, ascii_name(nm)); haskey(dsX, jn) || continue
            truth = r2X(jn); st = get!(istat, nm, Float64[0.0, 0.0])
            for k in 1:kmt
                lev = lev_of(k); fin(truth[c,lev]) || continue   # only skip cells MARBL itself leaves as FILL
                g = Gc[nm][1,1,k]
                # a non-finite Gⁿ where MARBL has a finite value is a FAILURE, never a silent skip
                st[2] += 1
                isfinite(g) || (st[1] = Inf; continue)
                # faithful-assembly check is interior-only: at the floor Gⁿ = pointwise transition + sediment/Δz
                k > 1 && (maxint_assembly = max(maxint_assembly, abs(g - biogeochemical_transition(1, 1, k, gc, bgcc, Val(nm), mc.clock, fc))))
                st[1] = max(st[1], abs(g - truth[c,lev]) / (abs(truth[c,lev]) + 1e-30))
            end
        end

        inj(nm) = Gc[nm][1,1,1] - biogeochemical_transition(1, 1, 1, gc, bgcc, Val(nm), mc.clock, fc)
        max_pocfloor = max(max_pocfloor, abs(sm.tracked_fields.POC_floor[1,1,1]))
        for nm in keys(min_inj); haskey(Gc, nm) && (min_inj[nm] = min(min_inj[nm], abs(inj(nm)))); end
        println("   $tag col $c/$ncX done (kmt=$kmt)"); flush(stdout)
    end
    close(dsX)

    @testset "$tag" begin
        @test maxint_assembly == 0.0                          # Gⁿ ≡ assembled pointwise biogeochemistry
        @test 0 < max_pocfloor < 1e-2                         # physical floor flux (flat grid gave 1e51 garbage)
        for nm in (:O₂,:NO₃,:PO₄,:Si,:DIC,:Alk); @test min_inj[nm] > 0; end   # sediment injects a dissolved return
        @test min_inj[:Fe] == 0                               # all particulate Fe is buried ⇒ no dissolved return
        for nm in sort(collect(keys(istat)); by = String)
            st = istat[nm]; st[2] == 0 && continue
            tol = nm === :Fe ? 1e-8 : 1e-9                    # J_Fe carries MARBL's cancellation-prone iron speciation
            @printf("    %-9s max-rel = %.2e over %4d cells\n", string(nm), st[1], Int(st[2]))
            @test st[1] < tol
        end
    end
end

@testset "stepped full model vs MARBL — Gⁿ (interior + sediment) at every cell, all 4 configs" begin
    stepped_config("base CESM (sp/diat/diaz)", HIST_base,  g -> MARBLPlankton(),         (:sp, :diat, :diaz))
    stepped_config("+cocco (4P1Z)",            HIST_cocco, g -> MARBLPlankton_cocco(g),  (:sp, :diat, :diaz, :cocco))
    stepped_config("general2 (2P1Z)",          HIST_gen("general2"), g -> general2_plankton(),   (:sp, :diat))
    stepped_config("general5 (5P2Z)",          HIST_gen("general5"), g -> general5_plankton(g),  (:sp, :diat, :diaz, :cocco, :diat2))
end

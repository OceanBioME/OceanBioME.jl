#####
##### Shared setup for the MARBL Fortran comparisons
#####
##### Where the baselines live, how to read them, and the pieces both harnesses need: the fill value the
##### driver writes below the sea floor, the tracer↔variable name map, and the carbon chemistry
##### configured to match MARBL's own solver.
#####

using NCDatasets
using Oceananigans.Fields: interior
using OceanBioME: CarbonChemistry, MorelMaritorenaPhotosyntheticallyActiveRadiation

const CCM = OceanBioME.Models.CarbonChemistryModel

#####
##### Where the data is
#####

const MARBL_ROOT   = normpath(joinpath(@__DIR__, "..", "..", "marbl-info3"))
const DRAFT_DIR    = joinpath(MARBL_ROOT, "draft")
const MARBL_INPUTS = joinpath(MARBL_ROOT, "MARBL", "tests", "input_files")

const HIST_BASE  = joinpath(MARBL_INPUTS, "baselines", "call_compute_subroutines.history.nc")
const HIST_COCCO = joinpath(MARBL_INPUTS, "baselines", "call_compute_subroutines.cocco.history.nc")
const INITIAL_CONDITIONS = joinpath(MARBL_INPUTS, "initial_conditions",
                                    "call_compute_subroutines.20190718.nc")

hist_general(name) = joinpath(DRAFT_DIR, "general_config", "baselines", "$name.history.nc")

"""
    baselines_available()

Whether the MARBL baseline files these comparisons read are present. They are several megabytes of
driver output which is not part of the package, so a harness reports and skips rather than failing
when they are absent.
"""
baselines_available() = all(isfile, (HIST_BASE, HIST_COCCO, INITIAL_CONDITIONS))

function require_baselines(name)
    baselines_available() && return true

    @warn "$name skipped: the MARBL baselines are not present. Expected them under $MARBL_INPUTS."

    return false
end

#####
##### Reading the driver's output
#####

# the driver writes this below the sea floor, so it marks a cell as outside the domain
const FILL = 9.96920996838687e36

fin(x) = (x != FILL) & isfinite(x)

median_of(v) = (s = sort(collect(v)); n = length(s); iseven(n) ? (s[n÷2] + s[n÷2+1]) / 2 : s[n÷2+1])

# a (column, level) array of one variable, with missing values turned into the fill value
read_column_variable(ds, name) = permutedims(Array{Float64}(coalesce.(ds[name][:, :], FILL)))

# MARBL's variable names for the tracers whose names differ from ours beyond the subscripts
const MARBL_VARIABLE = Dict(:NO₃ => "NO3", :NH₄ => "NH4", :PO₄ => "PO4", :Fe => "Fe", :Si => "SiO3",
                            :DIC => "DIC", :Alk => "ALK", :O₂ => "O2", :T => "insitu_temp")

marbl_variable(nm) = get(MARBL_VARIABLE, nm, ascii_name(nm))

#####
##### Carbon chemistry matched to MARBL's own solver
#####
##### Reproducing MARBL's carbonate system to machine precision needs three things together: its bracketed
##### safeguarded-Newton pH solve over pH ∈ [6, 9], the Dickson & Riley fluoride constant rather than
##### Perez & Fraga, and its gas constant. Any one of them alone leaves a residual.
#####

struct DicksonRileyKF{PC}
    pressure_correction :: PC
end

@inline (c::DicksonRileyKF)(T, S, Is, KS; P = nothing) =
    c.pressure_correction(T, P) *
    exp(1590.2 / T - 12.641 + 1.525 * sqrt(Is) + log(1 - 0.001005 * S) +
        log(1 + (0.14 / 96.062 * S / 1.80655) / KS))

const CC_MARBL = CarbonChemistry(;
    density_function = (args...) -> 1026.0,
    solver = OceanBioME.DRTSafeSolver(; max_iters = 100, bracket_grows = 3, atol = 1e-10,
                                      H_lower = 1e-9, H_upper = 1e-6),
    fluoride = DicksonRileyKF(CCM.PressureCorrection(; a₀ = -9.78, a₁ = -0.0090, a₂ = -0.000942,
                                                     b₀ = -0.00391, b₁ = 0.000054)),
    carbonic_acid = (
        K1 = CCM.K1(; pressure_correction = CCM.PressureCorrection(; a₀ = -25.50, a₁ = 0.1271, a₂ = 0.0,
                                                                   b₀ = -0.00308, b₁ = 0.0000877, R = 83.1451)),
        K2 = CCM.K2(; pressure_correction = CCM.PressureCorrection(; a₀ = -15.82, a₁ = -0.0219, a₂ = 0.0,
                                                                   b₀ = 0.00113, b₁ = -0.0001475, R = 83.1451))))

#####
##### Light
#####

# the chlorophyll dependent attenuation the driver uses, per cell
attenuation_per_cell(chl, Δz) =
    (w = max(chl, 0.02); (w < 0.13224 ? 0.000919 * w^0.3536 : 0.001131 * w^0.4562) * (Δz * 100))

# MARBL's f_qsw_par (marbl_settings_mod): the fixed fraction of surface shortwave that is PAR. The light
# model is now supplied the surface PAR directly, so this scaling (part of MARBL's surface sub-column
# handling) is applied here, in the harness, before the model sees it.
const F_QSW_PAR = 0.45

# MARBL's surface forcing is the shortwave irradiance in each of its ice radiation sub columns; a single
# column sees their area weighted mean, from which the light model makes the PAR itself
column_shortwave(areas, shortwave) =
    [sum(@view(areas[:, c]) .* @view(shortwave[:, c])) for c in axes(areas, 2)]

# each MARBL column is a separate i index on the comparison grid, so the surface forcing is read from a
# vector of one PAR per column
@inline column_surface_PAR(i, j, grid, clock, fields, Q) = @inbounds Q[i]

"""
    marbl_light(grid, shortwave)

The light model MARBL's `compute_PAR` implements, forced by one surface `shortwave` per column (scaled to
PAR by `F_QSW_PAR` here), which computes both the cell mean PAR and the interface PAR from the chlorophyll
(rather than the harnesses prescribing MARBL's `PAR_avg` and reconstructing an interface from it).
"""
marbl_light(grid, shortwave) =
    MorelMaritorenaPhotosyntheticallyActiveRadiation(grid, column_surface_PAR;
                                                     discrete_form = true, parameters = F_QSW_PAR .* shortwave)

# a computed light field back in MARBL's (column, level) top down indexing
par_top_down(field) = (A = Array(interior(field))[:, 1, :];
                       [A[c, size(A, 2) - lev + 1] for c in axes(A, 1), lev in axes(A, 2)])

"""
    par_deviation(field, PAR_avg, kmt; threshold = 1e-12)

How far a computed light `field` sits from MARBL's `PAR_avg`, as the largest absolute difference over
every water column cell and the largest relative difference over the cells holding more than
`threshold` W/m². The relative measure needs the floor because MARBL propagates each ice radiation sub
column separately and so cuts off at a different depth in each, which a single column cannot reproduce;
the cells that disagree hold ~1e-19 W/m², where a relative measure is meaningless.
"""
function par_deviation(field, PAR_avg, kmt; threshold = 1e-12)
    ours = par_top_down(field)

    difference = 0.0
    relative = 0.0

    for c in axes(ours, 1), lev in 1:kmt[c]
        fin(PAR_avg[c, lev]) || continue

        error = abs(ours[c, lev] - PAR_avg[c, lev])

        difference = max(difference, error)

        PAR_avg[c, lev] > threshold && (relative = max(relative, error / PAR_avg[c, lev]))
    end

    return (; difference, relative)
end

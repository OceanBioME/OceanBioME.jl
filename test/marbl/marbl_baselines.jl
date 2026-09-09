#####
##### Shared setup for the MARBL Fortran comparisons
#####
##### Where the baselines live, how to read them, and the pieces both harnesses need: the fill value the
##### driver writes below the sea floor, the tracer↔variable name map, and the carbon chemistry
##### configured to match MARBL's own solver.
#####

using NCDatasets
using OceanBioME: CarbonChemistry

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

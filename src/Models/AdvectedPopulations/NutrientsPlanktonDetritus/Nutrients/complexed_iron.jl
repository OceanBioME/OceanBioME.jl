using Oceananigans.Advection: advective_tracer_flux_z, Centered, ExplicitTimeDiscretization
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity
using Oceananigans.Operators: Δzᶜᶜᶜ, Azᶜᶜᶠ
using Oceananigans.Units: day

using OceanBioME.Light: @subcolumn_average, interface_par

struct ComplexedIron{FT}
          ligand_stability_constant :: FT  # L/µmol, free↔bound equilibrium; MARBL KFeLig1
               iron_scavenging_rate :: FT  # 1/yr per unit sinking mass; MARBL parm_Fe_scavenge_rate0
         bound_iron_scavenging_rate :: FT  # 1/yr per unit sinking mass; MARBL parm_FeLig_scavenge_rate0
             ligand_scavenging_rate :: FT  # 1/yr per unit sinking mass; MARBL parm_Lig_scavenge_rate0

            ligand_production_ratio :: FT  # ligand per unit remineralisation; MARBL remin_to_Lig
            ligand_degradation_rate :: FT  # bacterial degradation per unit POC remin; MARBL parm_Lig_degrade_rate0
        ligand_iron_uptake_fraction :: FT  # share of autotroph Fe uptake that removes ligand
    ligand_photochemistry_reference :: FT  # 1/(3 months) in 1/yr units; the surface photo-oxidation scale

    particulate_organic_matter_mass :: FT  # g/mol ballast mass per POC; MARBL 3.0·POC%mass
                       calcite_mass :: FT  # g/mol; MARBL P_CaCO3%mass
                          opal_mass :: FT  # g/mol; MARBL P_SiO2%mass
              dust_scavenging_scale :: FT  # mg/g for the dust ballast term
           sinking_mass_unit_scale  :: FT  # puts the SI sinking mass back into the units the rates are calibrated in

                 uv_reference_depth :: FT  # m, the depth the surface photo-oxidation dose is spread over
end

"""
    ComplexedIron([FT = Float64;] kwargs...)

The dissolved-iron system for the `iron` slot of [`Nutrients`](@ref), adding the tracers `Fe` (total
dissolved iron) and `Lig` (iron-binding ligand). It is a drop-in replacement for the plain `Fe`
nutrient: growth limitation is unchanged, because that depends on total dissolved iron.

Over plain `Fe` it adds single-ligand speciation into free and ligand-bound iron, scavenging of
dissolved iron onto sinking particles (a sink on `Fe` and a source into sinking particulate iron, so
the two stay conserved), and a prognostic ligand budget with production from remineralisation and loss
to scavenging, uptake, surface photo-oxidation and bacterial degradation.

Keyword Arguments
=================

- `ligand_stability_constant`: the free↔bound iron equilibrium constant (L/µmol)
- `iron_scavenging_rate`, `bound_iron_scavenging_rate`, `ligand_scavenging_rate`: scavenging base rates
  (1/yr per unit sinking mass)
- `ligand_production_ratio`, `ligand_degradation_rate`, `ligand_iron_uptake_fraction`: the ligand budget
- `particulate_organic_matter_mass`, `calcite_mass`, `opal_mass`, `dust_scavenging_scale`: the ballast
  masses that set how much sinking material is available to scavenge onto
- `uv_reference_depth`: the depth over which the surface photo-oxidation dose is spread (m)
"""
function ComplexedIron(FT = Float64;
                       ligand_stability_constant       = 10.0e13 * 1.0e-6,
                       iron_scavenging_rate            = 22.0,
                       bound_iron_scavenging_rate      = 1.2,
                       ligand_scavenging_rate          = 0.015,
                       ligand_production_ratio         = 1.0e-4,
                       ligand_degradation_rate         = 9.4e-5,
                       ligand_iron_uptake_fraction     = 0.20,
                       ligand_photochemistry_reference = 12 / 3,
                       particulate_organic_matter_mass = 3.0 * 12.01,   # g/mol
                       calcite_mass                    = 100.09,        # g/mol
                       opal_mass                       = 60.08,         # g/mol
                       dust_scavenging_scale           = 1.0e3,         # mg/g
                       sinking_mass_unit_scale         = 1.0e2,
                       uv_reference_depth              = 10.0)          # m

    return ComplexedIron{FT}(convert(FT, ligand_stability_constant),
                             convert(FT, iron_scavenging_rate),
                             convert(FT, bound_iron_scavenging_rate),
                             convert(FT, ligand_scavenging_rate),
                             convert(FT, ligand_production_ratio),
                             convert(FT, ligand_degradation_rate),
                             convert(FT, ligand_iron_uptake_fraction),
                             convert(FT, ligand_photochemistry_reference),
                             convert(FT, particulate_organic_matter_mass),
                             convert(FT, calcite_mass),
                             convert(FT, opal_mass),
                             convert(FT, dust_scavenging_scale),
                             convert(FT, sinking_mass_unit_scale),
                             convert(FT, uv_reference_depth))
end

required_biogeochemical_tracers(::ComplexedIron) = (:Fe, :Lig)
required_biogeochemical_auxiliary_fields(::ComplexedIron) = (:PAR, )

const ComplexedIron_NPD = NutrientsPlanktonDetritus{<:Any, <:Nutrients{<:Any, <:Any, <:ComplexedIron}}

# the scheme used to reconstruct the sinking ballast flux that drives scavenging; any consistent scheme
# infers the flux to its own order
const sinking_advection = Centered()

#####
##### Pure rate helpers
#####

# single-ligand speciation: total dissolved iron → free iron and ligand-bound iron, the positive root of
# 0 = K·Fe′² + (1 + K(L − Fe))·Fe′ − Fe. Non-positive iron gives zero of both; the square-root argument
# stays non-negative because the linear coefficient dominates.
@inline function iron_speciation(Fe, Lig, K)
    p₁ = one(Fe) + K * (Lig - Fe)

    Fefree = (-p₁ + sqrt(p₁^2 + 4K * Fe)) / (2K)
    Fefree = ifelse(Fe > zero(Fe), Fefree, zero(Fe))

    FeLig = K * Fefree * Lig / (one(Fe) + K * Fefree)

    return Fefree, FeLig
end

# the sinking mass available to scavenge onto, per unit area. The scavenging rate constants are calibrated
# against a mass built from fluxes in cgs, so `sinking_mass_unit_scale` restores those units and lets the
# rate constants be used unmodified; `dust_scavenging_scale` is a plain g→mg conversion for the same reason.
@inline sinking_mass(Fpoc, Fcaco3, Fsio2, Fdust, ci::ComplexedIron) =
    (Fpoc   * ci.particulate_organic_matter_mass +
     Fcaco3 * ci.calcite_mass +
     Fsio2  * ci.opal_mass +
     Fdust  * ci.dust_scavenging_scale) * ci.sinking_mass_unit_scale

# base rates are per year per unit sinking mass
@inline iron_scavenging_flux(Fefree, FeLig, M, ci::ComplexedIron) =
    (Fefree * ci.iron_scavenging_rate + FeLig * ci.bound_iron_scavenging_rate) * M / (365day)

@inline ligand_scavenging_flux(FeLig, M, ci::ComplexedIron) =
    FeLig * ci.ligand_scavenging_rate * M / (365day)

@inline function ligand_photochemistry_rate(i, j, k, grid, ci::ComplexedIron, aux)
    I       = interface_par(i, j, k + 1, grid, aux) # returns SubcolumnValues if available
    Δz      = Δzᶜᶜᶜ(i, j, k, grid)
    surface = k == grid.Nz

    return ligand_photochemistry_rate(I, ci.uv_reference_depth / Δz, ci.ligand_photochemistry_reference / (365day), surface)
end

@subcolumn_average @inline ligand_photochemistry_rate(PAR, uv_depth, rate, surface) =
    ifelse(surface & (PAR > 1),
           (log(max(PAR, one(PAR))) * oftype(PAR, 0.4373)) * uv_depth * rate,
           zero(PAR))

#####
##### Coupling functions
#####

# throws away negative sinking, maybe we should use the correct advection scheme somehow
@inline function downward_sinking_flux(i, j, k, grid, w, c)
    kf = k + 1

    F = advective_tracer_flux_z(i, j, kf, grid, sinking_advection, ExplicitTimeDiscretization(), w, c) /
        Azᶜᶜᶠ(i, j, kf, grid)

    return max(-F, zero(F))
end

# method for explicit particles
@inline function scavenging_sinking_mass(i, j, k, grid, bgc, fields, ci::ComplexedIron)
    wpoc  = biogeochemical_drift_velocity(bgc, Val(:POC)).w
    wbsi  = biogeochemical_drift_velocity(bgc, Val(:bSi)).w
    wcaco = biogeochemical_drift_velocity(bgc, Val(:CaCO₃)).w

    Fpoc  = downward_sinking_flux(i, j, k, grid, wpoc,  fields.POC)
    Fbsi  = downward_sinking_flux(i, j, k, grid, wbsi,  fields.bSi)
    Fcaco = downward_sinking_flux(i, j, k, grid, wcaco, fields.CaCO₃)

    return sinking_mass(Fpoc, Fcaco, Fbsi, zero(Fpoc), ci)
end

@inline function iron_scavenging(i, j, k, grid, ci::ComplexedIron, bgc, fields, aux)
    @inbounds begin
        Fe  = fields.Fe[i, j, k]
        Lig = fields.Lig[i, j, k]
    end

    Fefree, FeLig = iron_speciation(Fe, Lig, ci.ligand_stability_constant)

    M = scavenging_sinking_mass(i, j, k, grid, bgc, fields, ci)

    return iron_scavenging_flux(Fefree, FeLig, M, ci)
end

@inline ligand_production(i, j, k, grid, ci::ComplexedIron, bgc, fields, aux) =
    ci.ligand_production_ratio *
    (poc_remineralisation(i, j, k, grid, bgc.detritus, fields) +
     doc_production(i, j, k, grid, bgc, fields, aux))

@inline function ligand_scavenging(i, j, k, grid, ci::ComplexedIron, bgc, fields, aux)
    @inbounds begin
        Fe  = fields.Fe[i, j, k]
        Lig = fields.Lig[i, j, k]
    end

    _, FeLig = iron_speciation(Fe, Lig, ci.ligand_stability_constant)

    M = scavenging_sinking_mass(i, j, k, grid, bgc, fields, ci)

    return ligand_scavenging_flux(FeLig, M, ci)
end

@inline ligand_degradation(i, j, k, grid, ci::ComplexedIron, bgc, fields) =
    ci.ligand_degradation_rate * poc_remineralisation(i, j, k, grid, bgc.detritus, fields)

@inline ligand_photochemistry(i, j, k, grid, ci::ComplexedIron, bgc, fields, aux) =
    @inbounds fields.Lig[i, j, k] * ligand_photochemistry_rate(i, j, k, grid, ci, aux)

@inline ligand_loss(i, j, k, grid, ci::ComplexedIron, bgc, fields, aux) =
    ligand_scavenging(i, j, k, grid, ci, bgc, fields, aux) +
    ci.ligand_iron_uptake_fraction *
        nutrient_uptake(i, j, k, grid, Val(:Fe), bgc.plankton, bgc, fields, aux) +
    ligand_photochemistry(i, j, k, grid, ci, bgc, fields, aux) +
    ligand_degradation(i, j, k, grid, ci, bgc, fields)

# dissolved iron needs no override: it picks up the scavenging sink through the particulate-iron production
@inline (bgc::ComplexedIron_NPD)(i, j, k, grid, ::Val{:Lig}, clock, fields, aux) =
    ligand_production(i, j, k, grid, bgc.nutrients.iron, bgc, fields, aux) -
    ligand_loss(i, j, k, grid, bgc.nutrients.iron, bgc, fields, aux)

#####
##### Ballast coupling
#####
##### The ballast sweep already knows the sinking mass exactly, so scavenging reads it rather than
##### reconstructing it from sinking velocities that do not exist there, and the sweep's own mass flux is
##### this cycle's `sinking_mass`.
#####

# without an iron cycle nothing scavenges, so the ballast sweep's mass flux is unused
@inline ballast_sinking_mass(FPOC, FCaCO₃, FbSi, Fdust, iron, b) = zero(FPOC)

@inline ballast_sinking_mass(FPOC, FCaCO₃, FbSi, Fdust, ci::ComplexedIron, b) =
    sinking_mass(FPOC, FCaCO₃, FbSi, Fdust, ci)

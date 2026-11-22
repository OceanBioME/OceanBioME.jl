"""
    PrognosticLigand(FT = Float64; production_ratio = 1e-4,
                                   dissolved_production_ratio = production_ratio,
                                   degradation_ratio = 9.4e-5,
                                   uptake_ratio = 0.4,
                                   scavenging_rate = 0.015,
                                   photodegradation_rate = 0.4373 * 0.018 / 86400,
                                   photodegradation_depth = 10,
                                   low_ligand_threshold = 1.25e-3)

MARBL's prognostic iron-binding ligand (`compute_Lig_terms` and `compute_scavenging` in
`marbl_interior_tendency_mod.F90`) for the `total_ligand` of [`SimpleIron`](@ref), which adds the tracer `Lig`
(mmol m⁻³). The iron is partitioned between free and ligand-bound iron with the local `Lig`, and

    ∂Lig/∂t = production_ratio × Rₚ + dissolved_production_ratio × Dₚ
            - scavenging_rate × M × FeL
            - photodegradation_rate × ln(PAR) × (photodegradation_depth / Δz) × Lig     (top cell only, PAR > 1 W m⁻²)
            - min(1, Lig / low_ligand_threshold) × (uptake_ratio × U + degradation_ratio × Rₚ),

where `Rₚ` is the remineralisation of the sinking particulate organic matter and `Dₚ` the plankton's production of
dissolved organic matter (both in carbon, mmol C m⁻³ s⁻¹), `U` the plankton's gross iron uptake (mmol Fe m⁻³ s⁻¹; with
an `IronQuota` the quota times all of their production, as MARBL's `photoFe`, although the iron of the dissolved waste
is returned at once), `M` the
mass flux of the scavenging sinking particles (ng cm⁻² s⁻¹, as for the iron scavenging, so `SimpleIron` must scavenge
onto particles, e.g. with [`ParticleScavenging`](@ref)), and `PAR` (W m⁻²) that of the top cell.

MARBL's `DOC_prod` is only the semi-labile dissolved organic carbon: `parm_labile_ratio` (0.94) of the plankton's
dissolved losses goes straight to DIC, and the rest remineralises over 6–15 years. Where `Dₚ` includes labile matter
(e.g. a dissolved organic pool with a lifetime of months), `dissolved_production_ratio` of about 0.06 ×
`production_ratio` gives MARBL's ligand production per carbon produced.

Keyword Arguments
=================

- `production_ratio`: ligand produced per carbon of particulate organic matter remineralised, mol / mol C (MARBL's `remin_to_Lig`)
- `dissolved_production_ratio`: ligand produced per carbon made into dissolved organic matter, mol / mol C (MARBL's `remin_to_Lig` on its semi-labile `DOC_prod`)
- `degradation_ratio`: ligand degraded by bacteria per carbon remineralised, mol / mol C (`parm_Lig_degrade_rate0`)
- `uptake_ratio`: ligand taken up with the iron, mol / mol Fe
- `scavenging_rate`: of the ligand-bound iron's ligand, yr⁻¹ (ng cm⁻² s⁻¹)⁻¹ (`parm_Lig_scavenge_rate0`)
- `photodegradation_rate`: s⁻¹ per ln(W m⁻²)
- `photodegradation_depth`: the depth over which the top cell's photodegradation acts, m
- `low_ligand_threshold`: below which the uptake and bacterial degradation fall in proportion to the ligand, mmol m⁻³

```jldoctest
using OceanBioME

PrognosticLigand()

# output
PrognosticLigand{Float64} (:Lig)
├── Production: 0.0001 mol / mol C remineralised, 0.0001 mol / mol C dissolved
├── Bacterial degradation: 9.4e-5 mol / mol C remineralised
├── Uptake: 0.4 mol / mol Fe
├── Scavenging: 0.015 yr⁻¹ (ng cm⁻² s⁻¹)⁻¹ of ligand-bound iron
├── Photodegradation: 9.110416666666668e-8 s⁻¹ ln(PAR)⁻¹ over 10.0 m
└── Low ligand threshold: 0.00125 mmol/m³
```
"""
struct PrognosticLigand{FT}
              production_ratio :: FT # mol / mol C
    dissolved_production_ratio :: FT # mol / mol C
             degradation_ratio :: FT # mol / mol C
                  uptake_ratio :: FT # mol / mol Fe
               scavenging_rate :: FT # yr⁻¹ (ng cm⁻² s⁻¹)⁻¹
         photodegradation_rate :: FT # s⁻¹ ln(W m⁻²)⁻¹
        photodegradation_depth :: FT # m
          low_ligand_threshold :: FT # mmol m⁻³
end

PrognosticLigand(FT = Float64; production_ratio = 1e-4,
                               dissolved_production_ratio = production_ratio,
                               degradation_ratio = 9.4e-5,
                               uptake_ratio = 0.4,
                               scavenging_rate = 0.015,
                               photodegradation_rate = 0.4373 * 0.018 / 86400,
                               photodegradation_depth = 10,
                               low_ligand_threshold = 1.25e-3) =
    PrognosticLigand{FT}(production_ratio, dissolved_production_ratio, degradation_ratio, uptake_ratio, scavenging_rate,
                         photodegradation_rate, photodegradation_depth, low_ligand_threshold)

convert_ligand(FT, Lₜ::PrognosticLigand) =
    PrognosticLigand(FT; Lₜ.production_ratio, Lₜ.dissolved_production_ratio, Lₜ.degradation_ratio, Lₜ.uptake_ratio,
                         Lₜ.scavenging_rate, Lₜ.photodegradation_rate, Lₜ.photodegradation_depth, Lₜ.low_ligand_threshold)

ligand_summary(::PrognosticLigand) = "prognostic (:Lig)"
ligand_tracers(::PrognosticLigand) = (:Lig, )
@inline ligand_concentration(i, j, k, ::PrognosticLigand, fields) = @inbounds max(0, fields.Lig[i, j, k])

Base.summary(Lₜ::PrognosticLigand) = string(typeof(Lₜ), " (:Lig)")

function Base.show(io::IO, Lₜ::PrognosticLigand)
    print(io, summary(Lₜ), "\n",
          "├── Production: ", Lₜ.production_ratio, " mol / mol C remineralised, ",
                              Lₜ.dissolved_production_ratio, " mol / mol C dissolved\n",
          "├── Bacterial degradation: ", Lₜ.degradation_ratio, " mol / mol C remineralised\n",
          "├── Uptake: ", Lₜ.uptake_ratio, " mol / mol Fe\n",
          "├── Scavenging: ", Lₜ.scavenging_rate, " yr⁻¹ (ng cm⁻² s⁻¹)⁻¹ of ligand-bound iron\n",
          "├── Photodegradation: ", Lₜ.photodegradation_rate, " s⁻¹ ln(PAR)⁻¹ over ", Lₜ.photodegradation_depth, " m\n",
          "└── Low ligand threshold: ", Lₜ.low_ligand_threshold, " mmol/m³")
end

const PrognosticLigandNPD{FT} = NutrientsPlanktonDetritus{FT, <:Nutrients{<:Any, <:Any, <:SimpleIron{<:Any, <:Any, <:Any, <:PrognosticLigand}}}

@inline function (bgc::PrognosticLigandNPD{FT})(i, j, k, grid, ::Val{:Lig}, clock, fields, auxiliary_fields) where FT
    iron = bgc.nutrients.iron
    Lₜ = iron.total_ligand
    plankton = bgc.plankton

    R_CP = carbon_ratio(i, j, k, grid, plankton, bgc, fields)
    Rₚ = particulate_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields) * R_CP
    Dₚ = dissolved_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) * R_CP
    U  = gross_iron_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

    M = scavenging_particle_mass_flux(i, j, k, grid, iron.scavenging, bgc, fields)
    FeL = last(free_and_bound_iron(i, j, k, iron, fields))

    Lig = @inbounds fields.Lig[i, j, k]
    PAR = @inbounds auxiliary_fields.PAR[i, j, k]

    production = Lₜ.production_ratio * Rₚ + Lₜ.dissolved_production_ratio * Dₚ
    scavenging = Lₜ.scavenging_rate * M * FeL / convert(FT, 365 * 86400) # MARBL's year
    photodegradation = (k == size(grid, 3)) * Lₜ.photodegradation_rate * log(max(1, PAR)) *
                       Lₜ.photodegradation_depth / Δzᶜᶜᶜ(i, j, k, grid) * Lig
    biological = clamp(Lig / Lₜ.low_ligand_threshold, 0, 1) * (Lₜ.uptake_ratio * U + Lₜ.degradation_ratio * Rₚ)

    return production - scavenging - photodegradation - biological
end

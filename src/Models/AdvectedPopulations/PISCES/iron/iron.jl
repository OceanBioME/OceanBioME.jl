module Iron

export SimpleIron, IronTendencyArgs

using Oceananigans.Units
using Oceananigans.Grids: znode, Center

using OceanBioME.Models.PISCESModel: PISCES, anoxia_factor

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: 
    aggregation_of_colloidal_iron, degradation

using OceanBioME.Models.PISCESModel.ParticulateOrganicMatter: 
    iron_scavenging, iron_scavenging_rate, bacterial_iron_uptake,
    specific_degradation_rate, edible_flux_rate, edible_iron_flux_rate

using OceanBioME.Models.PISCESModel.Phytoplankton: uptake

using OceanBioME.Models.PISCESModel.Zooplankton: 
    non_assimilated_iron, upper_trophic_dissolved_iron,
    bacteria_concentration, bacteria_activity

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers
import OceanBioME.Models.PISCESModel: free_iron

include("simple_iron.jl")

@inline function ligand_concentration(iron::SimpleIron, DOC)
    Lₜᵐᵃˣ = iron.maximum_ligand_concentration
    Lₜ = iron.dissolved_ligand_ratio * DOC - Lₜᵐᵃˣ

    return max(Lₜᵐᵃˣ, Lₜ)
end

@inline function ligand_aggregation(iron::SimpleIron, Fe, DOC, T, scavenging_rate)
    Lₜ = ligand_concentration(iron, DOC)
    Fe′ = free_iron(iron, Fe, DOC, T)
    ΔFe = max(zero(Fe), Fe - Lₜ)

    return iron.excess_scavenging_enhancement * scavenging_rate * ΔFe * Fe′
end

@inline function free_iron(::SimpleIron, Fe, DOC, T)
    ligands = max(0.6, 0.09 * (DOC + 40) - 3)
    K = exp(16.27 - 1565.7 / max(T + 273.15, 5))
    Δ = 1 + K * ligands - K * Fe

    return (-Δ + √(Δ^2 + 4K * Fe)) / 2K
end

@inline function free_iron(iron::SimpleIron, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    DOC = @inbounds fields.DOC[i, j, k]
    Fe  = @inbounds fields.Fe[i, j, k]
    T   = @inbounds fields.T[i, j, k]

    return free_iron(iron, Fe, DOC, T)
end

end # module
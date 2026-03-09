module Iron

export SimpleIron

using Oceananigans.Units

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: 
    aggregation_of_colloidal_iron, degredation

using OceanBioME.Models.PISCESModel.ParticulateOrganicMatter: 
    iron_scavenging, iron_scavenging_rate, bacterial_iron_uptake

using OceanBioME.Models.PISCESModel.Phytoplankton: uptake

using OceanBioME.Models.PISCESModel.Zooplankton: 
    non_assimilated_iron, upper_trophic_dissolved_iron

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers
import OceanBioME.Models.PISCESModel: free_iron

include("simple_iron.jl")

@inline function ligand_concentration(iron::SimpleIron, DOC)
    Lₜᵐᵃˣ = iron.maximum_ligand_concentration
    Lₜ = iron.dissolved_ligand_ratio * DOC - Lₜᵐᵃˣ

    return max(Lₜᵐᵃˣ, Lₜ)
end

@inline function free_iron(::SimpleIron, Fe, DOC, T)
    # maybe some of these numbers should be parameters
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

@inline function ligand_aggregation(iron::SimpleIron, scavenging_rate, Fe, total_ligand_concentration, free_iron_concentration)
    λ̄ = iron.excess_scavenging_enhancement
    excess_iron = max(zero(Fe), Fe - total_ligand_concentration)

    return λ̄ * scavenging_rate * excess_iron * free_iron_concentration
end

@inline function iron_tendency(small_particles, grazing_waste, upper_trophic_waste,
                               consumption, ligand_aggregation_loss, colloidal_aggregation,
                               scavenging, bacterial_uptake)
    return (small_particles + grazing_waste + upper_trophic_waste
            - consumption - ligand_aggregation_loss - colloidal_aggregation - scavenging - bacterial_uptake)
end

end # module

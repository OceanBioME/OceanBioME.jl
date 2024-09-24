module Iron

export SimpleIron

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: 
    aggregation_of_colloidal_iron, bacterial_iron_uptake

using OceanBioME.Models.PISCESModel.ParticulateOrganicMatter: 
    iron_scavenging, iron_scavenging_rate, bacterial_iron_uptake

using OceanBioME.Models.PISCESModel.Phytoplankton: uptake

using OceanBioME.Models.PISCESModel.Zooplankton: 
    non_assimilated_iron, upper_trophic_iron_waste

import OceanBioME.Models.PISCESModel: free_iron

include("simple_iron.jl")

@inline function free_iron(::SimpleIron, i, j, k, grid, bgc, clock, fields)
    DOC = @inbounds fields.DOC[i, j, k]
    Fe  = @inbounds  fields.Fe[i, j, k]

    # maybe some of these numbers should be parameters
    ligands = max(0.6, 0.09 * (DOC + 40) - 3)
    K = exp(16.27 - 1565.7 / max(T + 273.15, 5))
    Δ = 1 + K * ligands - K * Fe

    return (-Δ + √(Δ^2 + 4K * Fe)) / 2K
end

end # module
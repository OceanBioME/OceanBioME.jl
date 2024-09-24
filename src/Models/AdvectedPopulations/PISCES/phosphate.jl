module Phosphates

export Phosphate

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: degredation

using OceanBioME.Models.PISCESModel.Phytoplankton: total_production

using OceanBioME.Models.PISCESModel.Zooplankton: inorganic_excretion, upper_trophic_respiration

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

struct Phosphate end

required_biogeochemical_tracers(::Phosphate) = tuple(:PO₄)

const PISCESPhosphate = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Phosphate}

@inline function (bgc::PISCESPhosphate)(i, j, k, grid, val_name::Val{:PO₄}, clock, fields, auxiliary_fields)
    θ = bgc.phosphate_redfield_ratio

    phytoplankton_uptake = total_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    grazing_waste = inorganic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    respiration_product = upper_trophic_respiration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    remineralisation = degredation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return θ * (grazing_waste + respiration_product + remineralisation - phytoplankton_uptake)
end

end # module
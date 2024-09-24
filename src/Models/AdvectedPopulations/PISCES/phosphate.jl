module Phosphates

export Phosphate

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: degredation

using OceanBioME.Models.PISCESModel.Phytoplankton: total_production

using OceanBioME.Models.PISCESModel.Zooplankton: inorganic_excretion, upper_trophic_respiration

struct Phosphate end

const PISCESPhosphate = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, Phosphate}

@inline function (bgc::PISCESPhosphate)(i, j, k, grid, val_name::Val{:PO₄}, clock, fields)
    θ = bgc.phosphate_redfield_ratio

    phytoplankton_uptake = total_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields)

    grazing_waste = inorganic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)

    respiration_product = upper_trophic_respiration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)

    remineralisation = degredation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields)

    return θ * (grazing_waste + respiration_product + remineralisation - phytoplankton_uptake)
end

end # module
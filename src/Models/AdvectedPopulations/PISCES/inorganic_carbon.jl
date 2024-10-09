module InorganicCarbons

export InorganicCarbon

using Oceananigans.Units

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: degredation

using OceanBioME.Models.PISCESModel.ParticulateOrganicMatter: 
    calcite_production, calcite_dissolution

using OceanBioME.Models.PISCESModel.Phytoplankton: total_production

using OceanBioME.Models.PISCESModel.Zooplankton:
    inorganic_excretion, upper_trophic_respiration

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

"""
    InorganicCarbon

Default parameterisation for `DIC`` and `Alk`alinity evolution. 
"""
struct InorganicCarbon end

required_biogeochemical_tracers(::InorganicCarbon) = (:DIC, :Alk)

const PISCESCarbon = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:InorganicCarbon}

@inline function (bgc::PISCESCarbon)(i, j, k, grid, ::Val{:DIC}, clock, fields, auxiliary_fields)
    zooplankton_respiration = inorganic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    upper_trophic = upper_trophic_respiration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    remineralisation = degredation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    calcite_diss = calcite_dissolution(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    calcite_prod = calcite_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    consumption = total_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return (zooplankton_respiration + upper_trophic + remineralisation 
            + calcite_diss - calcite_prod - consumption)
end

@inline function (bgc::PISCESCarbon)(i, j, k, grid, val_name::Val{:Alk}, clock, fields, auxiliary_fields)
    θ = bgc.nitrogen_redfield_ratio

    nitrate_production = bgc(i, j, k, grid, Val(:NO₃), clock, fields, auxiliary_fields)
    ammonia_production = bgc(i, j, k, grid, Val(:NH₄), clock, fields, auxiliary_fields)
    calcite_production = bgc(i, j, k, grid, Val(:CaCO₃), clock, fields, auxiliary_fields)

    # I think there are typos in Aumount 2015 but this is what it should be ( I think ???)
    return ammonia_production - nitrate_production - 2 * calcite_production
end

end # module
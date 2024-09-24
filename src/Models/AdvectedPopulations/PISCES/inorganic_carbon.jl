module InorganicCarbons

export InorganicCarbon

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: degredation

using OceanBioME.Models.PISCESModel.ParticulateOrganicMatter: 
    calcite_production, calcite_dissolution

using OceanBioME.Models.PISCESModel.Phytoplankton: total_production

using OceanBioME.Models.PISCESModel.Zooplankton:
    inorganic_excretion, upper_trophic_respiration

"""
    InorganicCarbon

Default parameterisation for `DIC`` and `Alk`alinity evolution. 
"""
struct InorganicCarbon end

const PISCESCarbon = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, InorganicCarbon}

@inline function (bgc::PISCESCarbon)(i, j, k, grid, val_name::Val{:DIC}, clock, fields)
    zooplankton_respiration = inorganic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)

    upper_trophic_respiration =  upper_trophic_respiration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)

    remineralisation = degredation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields)

    calcite_dissolution = calcite_dissolution(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields)

    calcite_production = calcite_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields)

    consumption = total_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields)

    return (zooplankton_respiration + upper_trophic_respiration + remineralisation 
            + calcite_dissolution - calcite_production
            - consumption)
end

@inline function (bgc::PISCESCarbon)(i, j, k, grid, val_name::Val{:Alk}, clock, fields)
    θ = bgc.nitrogen_redfield_ratio

    nitrate_production = bgc(i, j, k, grid, Val(:NO₃), clock, fields)
    ammonia_production = bgc(i, j, k, grid, Val(:NH₄), clock, fields)
    calcite_production = bgc(i, j, k, grid, Val(:CaCO₃), clock, fields)

    # I think there are typos in Aumount 2015 but this is what it should be ( I think ???)
    return ammonia_production - nitrate_production - 2 * calcite_production
end

end # module
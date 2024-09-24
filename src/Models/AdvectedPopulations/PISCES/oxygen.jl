module OxygenModels

export Oxygen

using Oceananigans.Units

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: 
    oxic_remineralisation, anoxic_remineralisation

using OceanBioME.Models.PISCESModel.Nitrogen: nitrification, nitrogen_fixation

using OceanBioME.Models.PISCESModel.Phytoplankton: uptake

using OceanBioME.Models.PISCESModel.Zooplankton:
    inorganic_excretion, upper_trophic_respiration

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

@kwdef struct Oxygen{FT}
    ratio_for_respiration :: FT = 133/122 # mol O₂ / mol C
   ratio_for_nitrification :: FT =  32/122 # mol O₂ / mol C
end

required_biogeochemical_tracers(::Oxygen) = tuple(:O₂)

const PISCESOxygen = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Oxygen}

@inline function (bgc::PISCESOxygen)(i, j, k, grid, val_name::Val{:O₂}, clock, fields, auxiliary_fields)
    θ_resp   = bgc.oxygen.ratio_for_respiration
    θ_nitrif = bgc.oxygen.ratio_for_nitrification

    zooplankton = θ_resp * inorganic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    upper_trophic = θ_resp * upper_trophic_respiration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    remineralisation = ((θ_resp + θ_nitrif) * oxic_remineralisation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
                        + θ_resp * anoxic_remineralisation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields))

    ammonia_photosynthesis = θ_resp * uptake(bgc.phytoplankton, Val(:NH₄), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    nitrate_photosynthesis = (θ_resp + θ_nitrif) * uptake(bgc.phytoplankton, Val(:NO₃), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    # I think (?) that we need the redfield raito here since θ_nitrif is per oxygen
    nitrif = θ_nitrif * nitrification(bgc.nitrogen, i, j, k, grid, bgc, clock, fields, auxiliary_fields) / bgc.nitrogen_redfield_ratio

    fixation = θ_nitrif * nitrogen_fixation(bgc.nitrogen, i, j, k, grid, bgc, clock, fields, auxiliary_fields) / bgc.nitrogen_redfield_ratio

    return (ammonia_photosynthesis + nitrate_photosynthesis + fixation
            - zooplankton - upper_trophic - nitrif)
end

end # module
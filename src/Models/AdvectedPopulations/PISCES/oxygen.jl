module OxygenModels

export Oxygen

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.DissolvedOrganicMatter: 
    oxic_remineralisation, anoxic_remineralisation

using OceanBioME.Models.PISCESModel.Nitrogen: nitrifcation, nitrogen_fixation

using OceanBioME.Models.PISCESModel.Phytoplankton: uptake

using OceanBioME.Models.PISCESModel.Zooplankton:
    inorganic_excretion, upper_trophic_respiration

@kwdef struct Oxygen{FT}
    ratio_for_respiration :: FT = 133/122 # mol O₂ / mol C
   ratio_for_nitrifcation :: FT = 32/122  # mol O₂ / mol C
end

const PISCESOxygen = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Any, Oxygen}

@inline function (bgc::PISCESOxygen)(i, j, k, grid, val_name::Val{:O₂}, clock, fields)
    θ_resp = oxy.ratio_for_respiration
    θ_nitrif  = oxy.ratio_for_nitrifcation

    zooplankton = θ_resp * inorganic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)

    upper_trophic = θ_resp * upper_trophic_respiration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)

    remineralisation = ((θ_resp + θ_nitrif) * oxic_remineralisation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields)
                        + θ_resp * anoxic_remineralisaiton(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields))

    ammonia_photosynthesis = θ_resp * uptake(bgc.phytoplankton, Val(:NH₄), i, j, k, grid, bgc, clock, fields)
    nitrate_photosynthesis = (θ_resp + θ_nitrif) * uptake(bgc.phytoplankton, Val(:NO₃), i, j, k, grid, bgc, clock, fields)

    # I think (?) that we need the redfield raito here since θ_nitrif is per oxygen
    nitrif = θ_nitrif * nitrification(bgc.nitrogen, i, j, k, grid, bgc, clock, fields) / bgc.nitrogen_redfield_ratio

    fixation = θ_nitrif * nitrogen_fixation(bgc.nitrogen, i, j, k, grid, bgc, clock, fields) / bgc.nitrogen_redfield_ratio

    return (ammonia_photosynthesis + nitrate_photosynthesis + fixation
            - zooplankton - upper_trophic - nitrif)
end

end # module
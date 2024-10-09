module Silicates

export Silicate

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.ParticulateOrganicMatter: 
    particulate_silicate_dissolution

using OceanBioME.Models.PISCESModel.Phytoplankton: silicate_uptake

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

struct Silicate end

required_biogeochemical_tracers(::Silicate) = tuple(:Si)

const PISCESSilicate = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, <:Silicate}

@inline function (bgc::PISCESSilicate)(i, j, k, grid, ::Val{:Si}, clock, fields, auxiliary_fields)
    consumption = silicate_uptake(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    dissolution = particulate_silicate_dissolution(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return dissolution - consumption
end

end # module
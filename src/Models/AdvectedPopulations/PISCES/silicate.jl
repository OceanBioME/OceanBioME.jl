module Silicates

export Silicate

using OceanBioME.Models.PISCESModel: PISCES

using OceanBioME.Models.PISCESModel.ParticulateOrganicMatter: 
    particulate_silicate_dissolution

using OceanBioME.Models.PISCESModel.Phytoplankton: silicate_uptake

struct Silicate end

const PISCESSilicate = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:Any, Silicate}

@inline function (bgc::PISCESSilicate)(i, j, k, grid, val_name::Val{:Si}, clock, fields)
    consumption = silicate_uptake(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields)

    dissolution = particulate_silicate_dissolution(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields)

    return dissolution - consumption
end

end # module
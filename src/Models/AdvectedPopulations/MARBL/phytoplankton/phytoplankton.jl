module Phytoplankton

export SomePhytoplanktonModel

using Oceananigans.Units

using OceanBioME.Models.MARBLModel: MARBL

# get the functions we need from the other groups
#using OceanBioME.Models.MARBLModel.Zooplankton: grazing

import Oceananigans.Biogeochemistry: required_biogeochemical_tracers

struct SomePhytoplanktonModel{FT}
    parameter1 :: FT
end

@inline function (bgc::MARBL{<:SomePhytoplanktonModel})(i, j, k, grid, val_name::Val{:P}, clock, fields, auxiliary_fields)
    P = @inbounds fields.P[i, j, k]

    growth = 1
    death = 0.1
    getting_grazed = 0.9

    return growth * P - death * P - getting_grazed * P^2
end

end # module
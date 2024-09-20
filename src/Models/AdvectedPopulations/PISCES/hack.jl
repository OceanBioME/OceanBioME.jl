using Oceananigans.Biogeochemistry: required_biogeochemical_tracers, required_biogeochemical_auxiliary_fields, extract_biogeochemical_fields
using Oceananigans.Grids: xnode, ynode, znode

import Oceananigans.Biogeochemistry: biogeochemical_transition

@inline function biogeochemical_transition(i, j, k, grid, bgc::AbstractContinuousFormBiogeochemistry, val_tracer_name, clock, fields)
    tracer_names_to_extract = required_biogeochemical_tracers(bgc)
    auxiliary_names_to_extract = required_biogeochemical_auxiliary_fields(bgc)

    tracer_fields_ijk = extract_biogeochemical_fields(i, j, k, grid, fields, tracer_names_to_extract)
    auxiliary_fields_ijk = extract_biogeochemical_fields(i, j, k, grid, fields, auxiliary_names_to_extract)

    x = xnode(i, j, k, grid, Center(), Center(), Center())
    y = ynode(i, j, k, grid, Center(), Center(), Center())
    z = znode(i, j, k, grid, Center(), Center(), Center())

    return bgc(val_tracer_name, x, y, z, clock.time, tracer_fields_ijk..., auxiliary_fields_ijk...)
end
using Oceananigans: defaults
using Oceananigans.Architectures: architecture
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.Models: fields
using Oceananigans.Utils: launch!

using OceanBioME.Models: CarbonChemistryModel

const gravitational_acceleration = defaults.gravitational_acceleration

function compute_calcite_saturation!(carbon_chemistry, calcite_saturation, model)
    grid = model.grid

    arch = architecture(grid)

    launch!(arch, grid, :xyz, _compute_calcite_saturation!, carbon_chemistry, calcite_saturation, grid, fields(model))

    fill_halo_regions!(calcite_saturation, model.grid, model.clock, fields(model))

    return nothing
end

@kernel function _compute_calcite_saturation!(carbon_chemistry, calcite_saturation, grid, model_fields)
    i, j, k = @index(Global, NTuple)

    T = @inbounds model_fields.T[i, j, k]
    S = @inbounds model_fields.S[i, j, k]
    DIC = @inbounds model_fields.DIC[i, j, k]
    Alk = @inbounds model_fields.Alk[i, j, k]
    silicate = @inbounds model_fields.Si[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())

    # very rough - don't think we should bother integrating the actual density
    # this is a type stability issue
    P = abs(z) * gravitational_acceleration * 1026 / 100000

    @inbounds calcite_saturation[i, j, k] = CarbonChemistryModel.calcite_saturation(carbon_chemistry; DIC, T, S, Alk, P, silicate)
end

using Oceananigans.Fields: Field, OneField

using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, immersed_cell

calculate_bottom_indices(grid) = OneField()

@kernel function find_bottom_cell(grid, bottom_indices)
    i, j = @index(Global, NTuple)

    Nz = size(grid, 3)

    k_bottom = 1

    while immersed_cell(i, j, k_bottom, grid.underlying_grid, grid.immersed_boundary) && k_bottom < Nz
        k_bottom += 1
    end

    @inbounds bottom_indices[i, j] = k_bottom
end

function calculate_bottom_indices(grid::ImmersedBoundaryGrid)
    arch = architecture(grid)
    indices = Field{Center, Center, Nothing}(grid, Int)

    launch!(arch, grid, :xy, find_bottom_cell, grid, indices)

    return indices
end

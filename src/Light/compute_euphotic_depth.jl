@kernel function _compute_euphotic_depth!(euphotic_depth, PAR, grid, cutoff)
    i, j = @index(Global, NTuple)

    surface_PAR = @inbounds PAR[i, j, grid.Nz]

    @inbounds euphotic_depth[i, j] = -Inf

    for k in grid.Nz-1:-1:1
        if (PAR[i, j, k] < surface_PAR * cutoff) && isinf(euphotic_depth[i, j])
            euphotic_depth[i, j] = znode(i, j, k, grid, Center(), Center(), Center())
        end
    end
end

function compute_euphotic_depth!(euphotic_depth, PAR, cutoff = 1/1000)
    grid = PAR.grid
    arch = architecture(grid)

    launch!(arch, grid, :xy, _compute_euphotic_depth!, euphotic_depth, PAR, grid, cutoff)

    fill_halo_regions!(euphotic_depth)

    return nothing
end
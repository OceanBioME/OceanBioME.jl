using Oceananigans.Fields: ConstantField, ZeroField

@kernel function _compute_euphotic_depth!(euphotic_depth, PAR, grid, cutoff)
    i, j = @index(Global, NTuple)

    surface_PAR = @inbounds (PAR[i, j, grid.Nz] + PAR[i, j, grid.Nz + 1])/2

    @inbounds euphotic_depth[i, j, 1] = -Inf

    for k in grid.Nz-1:-1:1
        PARₖ = @inbounds PAR[i, j, k]

        # BRANCHING!
        if (PARₖ <= surface_PAR * cutoff) && isinf(euphotic_depth[i, j])
            # interpolate to find depth
            PARₖ₊₁ = @inbounds PAR[i, j, k + 1]

            zₖ = znode(i, j, k, grid, Center(), Center(), Center())

            zₖ₊₁ = znode(i, j, k + 1, grid, Center(), Center(), Center())

            @inbounds euphotic_depth[i, j, 1] = zₖ + (log(surface_PAR * cutoff) - log(PARₖ)) * (zₖ - zₖ₊₁) / (log(PARₖ) - log(PARₖ₊₁))
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

# fallback for box models
compute_euphotic_depth!(::ConstantField, args...) = nothing
compute_euphotic_depth!(::ZeroField, args...) = nothing
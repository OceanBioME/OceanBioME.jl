using KernelAbstractions.Extras.LoopInfo: @unroll 
using Oceananigans.Operators: volume
using Oceananigans.Grids: AbstractGrid, Bounded, Periodic
using Oceananigans.Fields: fractional_indices
using Oceananigans.Architectures: arch_array

@inline get_node(::Bounded, i, N) = min(max(i, 1), N)
@inline get_node(::Periodic, i, N) = ifelse(i < 1, N , ifelse(i > N, 1, i))

# This won't work on the GPU and I can't see anyway around it, I think I need to come up with some much cleaverer way of doing it
@inline function get_nearest_nodes(x, y, z, grid, loc)
    i, j, k = fractional_indices(x, y, z, loc, grid)

    # Convert fractional indices to unit cell coordinates 0 <= (ξ, η, ζ) <=1
    # and integer indices (with 0-based indexing).
    ξ, i = modf(i)
    η, j = modf(j)
    ζ, k = modf(k)

    if (ξ, η, ζ) == (0, 0, 0) #particle on grid point special case
        return ((Int(i+1), Int(j+1), Int(k+1), 1), ), 1.0
    else
        nodes = repeat([(1, 1, 1, NaN)], 8)
        _normfactor = 0.0
        @unroll for n = 1:8
            # Move around cube corners getting node indices (0 or 1) and distances to them
            # Distance is d when the index is 0, or 1-d when it is 1
            a = 0 ^ (1 + (-1) ^ n)
            di = 0 ^ abs(1 - a) + ξ * (-1) ^ a

            b = 0 ^ (1 + (-1) ^ floor(n/2))
            dj = 0 ^ abs(1 - b) + η * (-1) ^ b

            c = 0 ^ (1 + (-1) ^ floor(n/4))
            dk = 0 ^ abs(1 - c) + ζ * (-1) ^ c

            @inbounds nodes[n] = (Int(i+1)+a, 
                                  Int(j+1) + b, 
                                  Int(k+1) + c, 
                                  sqrt(di^2 + dj^2 + dk^2))

            _normfactor += 1 / sqrt(di^2 + dj^2 + dk^2)
        end
        return nodes, 1/_normfactor
    end
end
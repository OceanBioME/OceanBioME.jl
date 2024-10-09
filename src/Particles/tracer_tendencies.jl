using Oceananigans.Grids: AbstractGrid, Bounded, Periodic

@inline get_node(::Bounded,  i, N) = min(max(i, 1), N)
@inline get_node(::Periodic, i, N) = ifelse(i < 1, N, ifelse(i > N, 1, i))

# TODO: update this with giant kelp stuff
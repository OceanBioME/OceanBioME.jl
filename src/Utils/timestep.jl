using Oceananigans.Advection: cell_advection_timescaleᶜᶜᶜ
using Oceananigans.Grids: Center, znode, minimum_zspacing

@inline function column_diffusion_timescale(model)
    grid = model.grid

      t = model.clock.time
     zₖ = znode(1, 1, k, grid, Center(), Center(), Center())
    Δzₖ = zspacing(1, 1, k, grid, Center(), Center(), Center())

    Δz2_ν = zeros(grid.Nz)

    @inbounds for k in 1:grid.Nz
        Δz2_ν[k] = Δzₖ^2 / model.closure.κ[1](0.0, 0.0, zₖ, t) # assumes all tracer closures are the same and x/y invariant
    end

    return minimum(Δz2_ν)
end

@inline column_advection_timescale(model) = minimum_zspacing(model.grid) / maximum_sinking_velocity(model.biogeochemistry)

@inline sinking_advection_timescale(model) = min(minimum_zspacing(model.grid) / maximum_sinking_velocity(model.biogeochemistry), cell_advection_timescaleᶜᶜᶜ(model))

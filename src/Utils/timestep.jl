using Oceananigans: Center
using Oceananigans.Advection: cell_advection_timescaleᶜᶜᶜ
using Oceananigans.Grids: ZRegRectilinearGrid, znodes

@inline Δz(::Center, k, grid::ZRegRectilinearGrid) = grid.Δzᵃᵃᶜ
@inline Δz(::Center, k, grid) = @inbounds grid.Δzᵃᵃᶜ[k]

@inline function column_diffusion_timescale(model)
    z = znodes(Center, model.grid)
    t = model.clock.time
    
    Δz2_ν = zeros(model.grid.Nz)

    @inbounds for k in 1:model.grid.Nz
        Δz2_ν[k] = Δz(Center(), k, model.grid) ^ 2 / model.closure.κ[1](0.0, 0.0, z[k], t) # assumes all tracer closures are the same and x/y invariant
    end

    return minimum(Δz2_ν)
end

@inline column_advection_timescale(model) = minimum(model.grid.Δzᵃᵃᶜ) / maximum_sinking_velocity(model.biogeochemistry)

@inline sinking_adveciton_timescale(model) = min(minimum(model.grid.Δzᵃᵃᶜ) / maximum_sinking_velocity(model.biogeochemistry), cell_advection_timescaleᶜᶜᶜ(model))
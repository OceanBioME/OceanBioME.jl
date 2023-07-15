module Morel

using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device
using Oceananigans.Utils: launch!

@kernel function _update_PAR(PAR¹, PAR², PAR³, Chl, zₑᵤ, grid, t, params) 
    i, j = @index(Global, NTuple)

    Nz = grid.Nz
    FT = eltype(grid)
    point_sixtytwo = FT(0.62)

    PAR⁰ = params.PAR⁰(t)

    # first point below surface
    @inbounds begin
        z = grid.zᵃᵃᶜ[Nz]

        k₁ = params.a.B + params.A.B*Chl[Nz] + params.b.B + params.B.B * Chl[Nz]^point_sixtytwo
        k₂ = params.a.G + params.A.G*Chl[Nz] + params.b.G + params.B.G * Chl[Nz]^point_sixtytwo
        k₃ = params.a.R + params.A.R*Chl[Nz] + params.b.R + params.B.R * Chl[Nz]^point_sixtytwo

        PAR¹[i, j, Nz] = PAR⁰ * params.PAR_frac.B * exp(k₁ * z)
        PAR²[i, j, Nz] = PAR⁰ * params.PAR_frac.G * exp(k₂ * z)
        PAR³[i, j, Nz] = PAR⁰ * params.PAR_frac.R * exp(k₃ * z)
    end

    # the rest of the points
    for k in grid.Nz-1:-1:1
        @inbounds begin
            z = grid.zᵃᵃᶜ[k]
            dz = grid.zᵃᵃᶜ[k+1] - z

            k₁ = params.a.B + params.A.B * Chl[k] + params.b.B + params.B.B*Chl[k]^point_sixtytwo
            k₂ = params.a.G + params.A.G * Chl[k] + params.b.G + params.B.G*Chl[k]^point_sixtytwo
            k₃ = params.a.R + params.A.R * Chl[k] + params.b.R + params.B.R*Chl[k]^point_sixtytwo

            PAR¹[i, j, k] = PAR¹[i, j, k+1] * exp(-k₁ * dz)
            PAR²[i, j, k] = PAR²[i, j, k+1] * exp(-k₂ * dz)
            PAR³[i, j, k] = PAR³[i, j, k+1] * exp(-k₃ * dz)

            if sum(PAR¹[i, j, k] + PAR²[i, j, k] + PAR³[i, j, k]) / PAR⁰ > FT(0.001)
                zₑᵤ[i, j] = - z #this and mixed layer depth are positive in PISCES
            end
        end
    end
end

"""
    update!(sim, params)

Function to integrate light attenuation using [Morel1988](@cite) model which should be called
in a callback as often as possible.

Requires three PAR bands to be supplied as auxiliary fields: `PAR¹`, `PAR²`, and `PAR³`
as well as `Chlᴾ` and `Chlᴰ` chlorophyll tracer fields.
"""
function update(sim, params)
    PAR¹, PAR², PAR³ = sim.model.auxiliary_fields.PAR¹, sim.model.auxiliary_fields.PAR², sim.model.auxiliary_fields.PAR³
    Chl = (sim.model.tracers.Chlᴾ .+ sim.model.tracers.Chlᴰ).*10^6

    zₑᵤ = sim.model.auxiliary_fields.zₑᵤ

    launch!(sim.model.architecture, sim.model.grid, :xy, _update_PAR,
            PAR¹, PAR², PAR³, Chl, zₑᵤ, sim.model.grid, sim.model.clock.time, params)
end

const defaults = (
    a = (B = 0.0145, G = 0.0754, R = 0.500),#m⁻¹
    A = (B = 0.044, G = 0.007, R = 0.007),#m²mg⁻¹
    b = (B = 0.0050, G = 0.00173, R = 0.0007),
    B = (B = 0.375, G = 0.290, R = 0.239),
    PAR_frac = (B = 1/3, G = 1/3, R = 1/3)
)

end #module

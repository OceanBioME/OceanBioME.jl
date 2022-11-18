module twoBands
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device
using Oceananigans.Utils: launch!

@kernel function _update_PAR!(PAR, grid, P, t, params) 
    i, j = @index(Global, NTuple) 

    sp = params.surface_PAR(t)

    z = grid.zᵃᵃᶜ[grid.Nz]

    ∫chlᵉʳ = -(P[i, j, grid.Nz]*params.Rd_chl/params.r_pig)^params.e_r*z
    ∫chlᵉᵇ = -(P[i, j, grid.Nz]*params.Rd_chl/params.r_pig)^params.e_b*z

    # TODO: Change this to an interpolation to the midpoint between the cells/make this actually valid for all grid types
    PAR[i, j, grid.Nz] =  sp*(exp(params.k_r0 * z - params.Χ_rp * ∫chlᵉʳ).+ exp(params.k_b0 * z - params.Χ_bp * ∫chlᵉᵇ))/2
    for k=grid.Nz-1:-1:1
        z = grid.zᵃᵃᶜ[k]
        dz = grid.zᵃᵃᶜ[k+1] - z 

        mean_pig = (P[i, j, k]+P[i, j, k+1])*params.Rd_chl/(2*params.r_pig)
        ∫chlᵉʳ += (mean_pig^params.e_r)*dz
        ∫chlᵉᵇ += (mean_pig^params.e_b)*dz

        PAR[i, j, k] =  sp*(exp(params.k_r0 * z - params.Χ_rp * ∫chlᵉʳ) + exp(params.k_b0 * z - params.Χ_bp * ∫chlᵉᵇ))/2
    end
end 

"""
    Light.twoBands.update!(sim, params)

Function to integrate light attenuation using [Karleskind2011](@cite) model which should be called in a callback as often as possible.

Requires a single `PAR` field to be defined in the auxiliary fields and a single `P` phytoplankton class.
"""
function  update!(sim, params)
    par_calculation =  launch!(sim.model.architecture, sim.model.grid, :xy, _update_PAR!,
                                   sim.model.auxiliary_fields.PAR, sim.model.grid, sim.model.tracers.P, time(sim), params,
                                   dependencies = Event(device(sim.model.architecture)))

    wait(device(sim.model.architecture), par_calculation)
end

defaults = (
    k_r0 = 0.225,  # m⁻¹
    k_b0 = 0.0232,  # m⁻¹
    Χ_rp = 0.037,  # m⁻¹(mgChlm⁻³)⁻ᵉʳ
    Χ_bp = 0.074,  # m⁻¹(mgChlm⁻³)⁻ᵉᵇ
    e_r = 0.629, 
    e_b = 0.674, 
    r_pig = 0.7,
)

end
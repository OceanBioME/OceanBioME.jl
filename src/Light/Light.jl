#This can not be computed as a computed field since then it would have to integrate on every grid point
#A better solution than callbacks may exist
"
Light attenuation by chlorophyll as described by Karleskind et al. 2011 (implimented as 2λ)
- intend to change to model by Morel 2001 in the future

These should be used by setting up a callback like:
`simulation.callbacks[:update_par] = Callback(Light.update_2λ!, IterationInterval(1), merge(params, (surface_PAR=surface_PAR,)))`

References:
Karleskind, P., Lévy, M., and Memery, L. (2011), Subduction of carbon, nitrogen, and oxygen in the northeast Atlantic, J. Geophys. Res., 116, C02025, doi:10.1029/2010JC006446. 
Morel, A., and Maritorena, S. (2001). Bio-optical properties of oceanic waters: a reappraisal. J. Geophys. Res. Oceans 106, 7163–7180. doi: 10.1029/2000JC000319
"
module Light
using Oceananigans
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device

include("2band.jl")

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
#This can not be computed as a computed field since then it would have to integrate on every grid point
#A better solution than callbacks may exist

module Light
using Oceananigans
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device

@kernel function _update_2λ!(par, grid, chl, t, params) 
    i, j = @index(Global, NTuple) 
    zpar=Oceananigans.Operators.znodes(par)

    ∫chlᵉʳ = Oceananigans.Architectures.arch_array(grid.architecture, zeros(grid.Nz))#I've checked and this array conversion doesn't cost anything (on CPU at least)
    @inbounds ∫chlᵉʳ[grid.Nz] = (chl[i, j, grid.Nz]/params.r_pig)^params.e_r*zpar[grid.Nz]

    ∫chlᵉᵇ = Oceananigans.Architectures.arch_array(grid.architecture, zeros(grid.Nz))
    @inbounds ∫chlᵉᵇ[grid.Nz] = (chl[i, j, grid.Nz]/params.r_pig)^params.e_b*zpar[grid.Nz]

    for k=grid.Nz-1:-1:1
        mean_pig::Float64 = (chl[i, j, k]+chl[i, j, k+1])/(2*params.r_pig)
        @inbounds ∫chlᵉʳ[k] = ∫chlᵉʳ[k+1] + (mean_pig^params.e_r)*(zpar[k]-zpar[k+1])
        @inbounds ∫chlᵉᵇ[k] = ∫chlᵉᵇ[k+1] + (mean_pig^params.e_b)*(zpar[k]-zpar[k+1])
    end

    @inbounds par[i, j, 1:end-3] = params.surface_PAR(t) .* (exp.(params.k_r0 .* zpar .+ params.Χ_rp .* ∫chlᵉʳ) .+ exp.(params.k_b0 .* zpar .+ params.Χ_bp .* ∫chlᵉᵇ))./2
end 

function  update_2λ!(sim, params)
    par_calculation = Oceananigans.Utils.launch!(sim.model.architecture, sim.model.grid, :xy, _update_2λ!,
                                   sim.model.auxiliary_fields.PAR, sim.model.grid, sim.model.tracers.P .* params.Rd_chl, time(sim), params,
                                   dependencies = Event(device(sim.model.architecture)))

    wait(device(sim.model.architecture), par_calculation)
end
end
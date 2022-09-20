@kernel function _update_2λ!(PAR, grid, P, t, params) 
    i, j = @index(Global, NTuple) 

    sp = params.surface_PAR(t)

    z = grid.zᵃᵃᶜ[grid.Nz]

    ∫chlᵉʳ = -(P[i, j, grid.Nz]*params.Rd_chl/params.r_pig)^params.e_r*z
    ∫chlᵉᵇ = -(P[i, j, grid.Nz]*params.Rd_chl/params.r_pig)^params.e_b*z
    PAR[i, j, grid.Nz] =  sp*(exp(params.k_r0 * z - params.Χ_rp * ∫chlᵉʳ).+ exp(params.k_b0 * z - params.Χ_bp * ∫chlᵉᵇ))/2
    for k=grid.Nz-1:-1:1
        z = grid.zᵃᵃᶜ[k]
        dz = grid.zᵃᵃᶜ[k+1] - z 

        mean_pig = (P[i, j, k]+P[i, j, k+1])*params.Rd_chl/(2*params.r_pig)
        ∫chlᵉʳ += (mean_pig^params.e_r)*dz
        ∫chlᵉᵇ += (mean_pig^params.e_b)*dz

        PAR[i, j, k] =  sp*(exp(params.k_r0 * dz - params.Χ_rp * ∫chlᵉʳ) + exp(params.k_b0 * dz - params.Χ_bp * ∫chlᵉᵇ))/2
    end
end 

function  update_2λ!(sim, params)
    par_calculation = Oceananigans.Utils.launch!(sim.model.architecture, sim.model.grid, :xy, _update_2λ!,
                                   sim.model.auxiliary_fields.PAR, sim.model.grid, sim.model.tracers.P, time(sim), params,
                                   dependencies = Event(device(sim.model.architecture)))

    wait(device(sim.model.architecture), par_calculation)
end
@kernel function _update_2λ!(par, grid, P, t, params) 
    i, j = @index(Global, NTuple) 

    surface_PAR = params.surface_PAR(t)

    ∫chlᵉʳ = (@inbounds P[i, j, grid.Nz]*params.Rd_chl/params.r_pig)^params.e_r*znode(Center(), Center(), Center(), i, j, grid.Nz)
    ∫chlᵉᵇ = (@inbounds P[i, j, grid.Nz]*params.Rd_chl/params.r_pig)^params.e_b*znode(Center(), Center(), Center(), i, j, grid.Nz)
    @inbounds PAR[i, j, grid.Nz] =  surface_PAR*(exp(params.k_r0 * zpar + params.Χ_rp * ∫chlᵉʳ).+ exp(params.k_b0 * zpar + params.Χ_bp * ∫chlᵉᵇ))/2
    for k=grid.Nz-1:-1:1
        z = znode(Center(), Center(), Center(), i, j, k)
        dz = z - znode(Center(), Center(), Center(), i, j, k+1)

        mean_pig::Float64 = (@inbounds P[i, j, k]*params.Rd_chl+@inbounds P[i, j, k+1]*params.Rd_chl)/(2*params.r_pig)
        ∫chlᵉʳ += (mean_pig^params.e_r)*dz
        ∫chlᵉᵇ += (mean_pig^params.e_b)*dz

        @inbounds PAR[i, j, k] =  surface_PAR*(exp(params.k_r0 * z + params.Χ_rp * ∫chlᵉʳ).+ exp(params.k_b0 * z + params.Χ_bp * ∫chlᵉᵇ))/2
    end
end 

function  update_2λ!(sim, params)
    par_calculation = Oceananigans.Utils.launch!(sim.model.architecture, sim.model.grid, :xy, _update_2λ!,
                                   sim.model.auxiliary_fields.PAR, sim.model.grid, sim.model.tracers.P, time(sim), params,
                                   dependencies = Event(device(sim.model.architecture)))

    wait(device(sim.model.architecture), par_calculation)
end

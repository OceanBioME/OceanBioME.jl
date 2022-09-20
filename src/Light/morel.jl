module Morel
#TODO: Add diurnal variability
@kernel function _update(PAR₁, PAR₂, PAR₃, Chl, zₑᵤ, grid, t, params) 
    i, j = @index(Global, NTuple) 
    Nz = grid.Nz
            
    PAR₀ = params.PAR₀(t)
    z = grid.zᵃᵃᶜ[Nz]

    k₁ = params.a.B + params.A.B*Chl[Nz] + params.b.B + params.B.B*Chl[Nz]^0.62
    k₂ = params.a.G + params.A.G*Chl[Nz] + params.b.G + params.B.G*Chl[Nz]^0.62
    k₃ = params.a.R + params.A.R*Chl[Nz] + params.b.R + params.B.R*Chl[Nz]^0.62

    PAR₁[i, j, Nz] = PAR₀*params.PAR_frac.B*exp(k₁*z)
    PAR₂[i, j, Nz] = PAR₀*params.PAR_frac.G*exp(k₂*z)
    PAR₃[i, j, Nz] = PAR₀*params.PAR_frac.R*exp(k₃*z)

    for k=grid.Nz-1:-1:1
        z = grid.zᵃᵃᶜ[k]
        dz = grid.zᵃᵃᶜ[k+1] - z 

        k₁ = params.a.B + params.A.B*Chl[k] + params.b.B + params.B.B*Chl[k]^0.62
        k₂ = params.a.G + params.A.G*Chl[k] + params.b.G + params.B.G*Chl[k]^0.62
        k₃ = params.a.R + params.A.R*Chl[k] + params.b.R + params.B.R*Chl[k]^0.62
    
        PAR₁[i, j, k] = PAR₁[i, j, k+1]*exp(-k₁*dz)
        PAR₂[i, j, k] = PAR₂[i, j, k+1]*exp(-k₂*dz)
        PAR₃[i, j, k] = PAR₃[i, j, k+1]*exp(-k₃*dz)

        if sum(PAR₁[i, j, k] + PAR₂[i, j, k] + PAR₃[i, j, k])/PAR₀ > 0.001
            zₑᵤ[i, j] = - z #this and mixed layer depth are positive in PISCES
        end
    end
end 

function  update(sim, params)
    PAR₁, PAR₂, PAR₃ = sim.model.auxiliary_fields.PAR₁, sim.model.auxiliary_fields.PAR₂, sim.model.auxiliary_fields.PAR₃
    Chl = (sim.model.tracers.Chlᴾ .+ sim.model.tracers.Chlᴰ).*10^6

    zₑᵤ = sim.model.auxiliary_fields.zₑᵤ

    par_calculation = Oceananigans.Utils.launch!(sim.model.architecture, sim.model.grid, :xy, _update,
                                   PAR₁, PAR₂, PAR₃, Chl, zₑᵤ, sim.model.grid, sim.model.clock.time
                                   dependencies = Event(device(sim.model.architecture)))

    wait(device(sim.model.architecture), par_calculation)
end

const defaults = (
    a = (B = 0.0145, G = 0.0754, R = 0.500),#m⁻¹
    A = (B = 0.044, G = 0.007, R = 0.007),#m²mg⁻¹
    b = (B = 0.0050, G = 0.00173, R = 0.0007),
    B = (B = 0.375, G = 0.290, R = 0.239)
)
end
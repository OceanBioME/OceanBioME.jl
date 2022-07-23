module NPZ
# NPZ model based on: 
#The reduction of plankton biomass induced by mesoscale stirring: A modeling study in the Benguela upwelling
# Ismael Hernández-Carrascoa et al 
using Oceananigans

import BGC: BGCModel

#source functions
P_forcing(x, y, z, t, P, N, Z, params) = params.β*N*P/(params.κₙ+N)*exp(z/params.hₗ)-params.α*params.η*P^2*Z/(params.α+params.η*P^2) - params.μₚ*P
N_forcing(x, y, z, t, P, N, Z, params) = params.Sₙ*(params.Nᵦ-N) - params.β*N*P/(params.κₙ+N) + params.μₙ*((1-params.γ)*params.α*params.η*P^2*Z/(params.α+params.η*P^2)+params.μₚ*P + params.μZ*Z^2)
Z_forcing(x, y, z, t, P, N, Z, params) = params.γ*params.α*params.η*P^2*Z/(params.α+params.η*P^2) - params.μZ*Z^2

# Define parameters
NPZ_parameters = (β=1/day, # Phytoplankton growth rate
η=1/day, # Prey capture rate
γ=0.75, # Assimilation efficiency of zooplankton
α=2/day, # Maximum grazing rate
κₙ=0.55, # Half-saturation constant for N uptake
μₙ=0.2, # Inefficiency of remineralization
μₚ=0.03/day, # Specific mortality rate
μZ=0.2/day, # Zooplankton mortality
hₗ=5, # Light (growth) e-folding depth
Sₙ=1/day, # rate of nutrient restoring
Nᵦ=8) # Nutrient concentration below mixed layer

function setup(grid, parameters, forcings=(T=nothing, S=nothing, PAR=nothing), fluxboundaries=(dic=AirSeaFlux.dic, ))
    if keys(forcings) != (:T, :S, :PAR)
        throw(ArgumentError("Need to provide all external forcings (T, S, PAR), either pass 'nothing' for a field to be used, a function of x, y, z, t, or an externally defined field which must also be an auxililiary field of the model. Currently T and S must be the same types."))
    end

    #setup forcings

    function_parameters = parameters   #this can be an auxiliary field or a function of x, y, z, t

    N_RHS = Forcing(N_forcing, field_dependencies = (:N,:P,:Z), parameters = function_parameters)
    P_RHS = Forcing(P_forcing, field_dependencies = (:N,:P,:Z), parameters = function_parameters)
    Z_RHS = Forcing(Z_forcing, field_dependencies = (:N,:P,:Z), parameters = function_parameters)

    #setup boundary conditions
    N_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))
    P_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))
    Z_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))  

   return BGCModel((:N, :P, :Z, #tracers
            (N=N_RHS, P=P_RHS, Z=Z_RHS), #forcing
            (N=N_bcs, P=P_bcs, Z=Z_bcs)) #boundaries
end
end # module

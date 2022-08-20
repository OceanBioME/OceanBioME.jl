"
NPZ biogeochemistry model as described by Hernández-Carrasco et al. 2014 and earlier publications.

WARNING: this model currently is not correctly configured

Variables: N, P, Z

Forcing: PAR

References:
Hernández-Carrasco, I., Rossi, V., Hernández-García, E., Garçon, V. and López, C., 2014. The reduction of plankton biomass induced by mesoscale stirring: A modeling study in the Benguela upwelling. Deep Sea Research Part I: Oceanographic Research Papers, 83, pp.65-80.
"
module NPZ
# NPZ model based on: 
#The reduction of plankton biomass induced by mesoscale stirring: A modeling study in the Benguela upwelling
# Ismael Hernández-Carrascoa et al 
using Oceananigans
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

#source functions
P_forcing(x, y, z, t, P, N, Z, params) = params.β*N*P/(params.κₙ+N)*exp(z/params.hₗ)-params.α*params.η*P^2*Z/(params.α+params.η*P^2) - params.μₚ*P
N_forcing(x, y, z, t, P, N, Z, params) = params.Sₙ*(params.Nᵦ-N) - params.β*N*P/(params.κₙ+N) + params.μₙ*((1-params.γ)*params.α*params.η*P^2*Z/(params.α+params.η*P^2)+params.μₚ*P + params.μZ*Z^2)
Z_forcing(x, y, z, t, P, N, Z, params) = params.γ*params.α*params.η*P^2*Z/(params.α+params.η*P^2) - params.μZ*Z^2

# Define parameters
defaults = (β=1/day, # Phytoplankton growth rate
η=1/day, # Prey capture rate
γ=0.75, # Assimilation efficiency of zooplankton
α=2/day, # Maximum grazing rate
κₙ=0.55, # Half-saturation constant for N uptake
μₙ=0.2, # Inefficiency of remineralization
μₚ=0.03/day, # Specific mortality rate
μZ=0.2/day, # Zooplankton mortality
hₗ=5, # Light (growth) e-folding depth
Sₙ=1/day, # rate of nutrient restoring
Nᵦ=8,# Nutrient concentration below mixed layer
Rd_chl = 1.31) 

tracers=(:N, :P, :Z)
optional_tracers=NamedTuple()

forcing_functions=(N=N_forcing, P=P_forcing, Z=Z_forcing)
sinking=NamedTuple()
end # module

using AIBECS, OceanBioME
using Oceananigans.Units
using Oceananigans.Biogeochemistry: biogeochemical_drift_velocity

import OceanBioME: setup_velocity_fields

@inline setup_velocity_fields(drift_speeds, grid, open_bottom) = drift_speeds

# not sure what this does
import AIBECS: @initial_value, initial_value

grid, circulaiton = OCIM2.load()

N_transport(p), P_transport(p), Z_transport(p) = circulaiton, circulaiton, circulaiton
D_transport(p) = circulaiton#transportoperator(grid, z -> w(z,p))

#@inline w(z, p) = 2.7489/day # @inbounds biogeochemical_drift_velocity(p.biogeochemistry, :D) # change for non-constant z

@inline numerical_values(unitful_value) = unitful_value.val

function update_PAR(P, grid, params)
    return (params.PAR⁰ * exp.(-params.kʷ .* numerical_values.(grid.depth_3D)))[iswet(grid)]
end

function N_change(N, P, Z, D, params)
    bgc = params.biogeochemistry

    params.light_field .= update_PAR(P, params.grid, params.light_model)
    T = 12.0

    return bgc.(Val(:N), nothing, nothing, nothing, nothing, N, P, Z, D, T, params.light_field)
end

function P_change(N, P, Z, D, params)
    bgc = params.biogeochemistry

    T = 12.0

    return bgc.(Val(:P), nothing, nothing, nothing, nothing, N, P, Z, D, T, params.light_field)
end

function Z_change(N, P, Z, D, params)
    bgc = params.biogeochemistry

    T = 12.0

    return bgc.(Val(:Z), nothing, nothing, nothing, nothing, N, P, Z, D, T, params.light_field)
end

function D_change(N, P, Z, D, params)
    bgc = params.biogeochemistry

    T = 12.0

    return bgc.(Val(:D), nothing, nothing, nothing, nothing, N, P, Z, D, T, params.light_field)
end

@initial_value struct AdaptModel{BGC, PAR, G, F, FT} <: AbstractParameters{FT}
    biogeochemistry::BGC |  nothing
    light_model::PAR     |  (PAR⁰ = 100.0, kʷ = 0.1, field = zeros(1))
    light_field::F       |  1.0
    grid::G              |  1.0
    anum::FT             |  1.0
end

num_of_wet_boxes = sum(iswet(grid))

biogeochemistry = NutrientPhytoplanktonZooplanktonDetritus(; grid = nothing, light_attenuation_model = nothing)
light_field = zeros(num_of_wet_boxes)

model = AdaptModel(; biogeochemistry, grid, light_field)

F = AIBECSFunction((N_transport, P_transport, Z_transport, D_transport), (N_change, P_change, Z_change, D_change), num_of_wet_boxes)

P = zeros(4num_of_wet_boxes)

P[1:num_of_wet_boxes] .= 10.0
P[num_of_wet_boxes + 1:2num_of_wet_boxes] .= 1.0
P[2num_of_wet_boxes + 1:3num_of_wet_boxes] .= 0.1

problem = SteadyStateProblem(F, P, model)

#solution = solve(problem, CTKAlg()).u

#plothorizontalslice(P, grid, depth=0m, color=:algae, rev=true)
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device, architecture
using Oceananigans.Utils: launch!

import Oceananigans.Fields: compute!

@kernel function update_TwoBandPhotosyntheticallyActiveRatiation!(PAR, grid, t) 
    i, j = @index(Global, NTuple)

    data = PAR.data
    P = PAR.operand.P_field
    
    PAR⁰ = PAR.operand.surface_intensity(t)

    kʳ = PAR.operand.water_red_attenuation
    kᵇ = PAR.operand.water_blue_attenuation
    χʳ = PAR.operand.chlorophyll_red_attenuation
    χᵇ = PAR.operand.chlorophyll_blue_attenuation
    eʳ = PAR.operand.chlorophyll_red_exponent
    eᵇ = PAR.operand.chlorophyll_blue_exponent
    r = PAR.operand.pigment_ratio
    Rᶜₚ = PAR.operand.phytoplankton_chlorophyll_ratio

    zᶜ = znodes(Center, grid)
    zᶠ = znodes(Face, grid)
    
    ∫chlʳ = @inbounds - (zᶜ[grid.Nz] - zᶠ[grid.Nz]) * (P[i, j, grid.Nz] * Rᶜₚ / r) ^ eʳ
    ∫chlᵇ = @inbounds - (zᶜ[grid.Nz] - zᶠ[grid.Nz]) * (P[i, j, grid.Nz] * Rᶜₚ / r) ^ eᵇ
    @inbounds data[i, j, grid.Nz] =  PAR⁰ * (exp(kʳ * zᶜ[grid.Nz] - χʳ * ∫chlʳ) + exp(kᵇ * zᶜ[grid.Nz] - χᵇ * ∫chlᵇ)) / 2

    @inbounds for k in grid.Nz-1:-1:1
        ∫chlʳ += (zᶜ[k + 1] - zᶠ[k]) * (P[i, j, k+1] * Rᶜₚ / r) ^ eʳ + (zᶠ[k] - zᶜ[k]) * (P[i, j, k] * Rᶜₚ / r) ^ eʳ
        ∫chlᵇ += (zᶜ[k + 1] - zᶠ[k]) * (P[i, j, k+1] * Rᶜₚ / r) ^ eᵇ + (zᶠ[k] - zᶜ[k]) * (P[i, j, k] * Rᶜₚ / r) ^ eᵇ
        data[i, j, k] =  PAR⁰ * (exp(kʳ * zᶜ[k] - χʳ * ∫chlʳ) + exp(kᵇ * zᶜ[k] - χᵇ * ∫chlᵇ)) / 2
    end
end 

struct TwoBandPhotosyntheticallyActiveRatiationOperand{FT, SI, P, C}
    water_red_attenuation :: FT
    water_blue_attenuation :: FT
    chlorophyll_red_attenuation :: FT
    chlorophyll_blue_attenuation :: FT
    chlorophyll_red_exponent :: FT
    chlorophyll_blue_exponent :: FT
    pigment_ratio :: FT

    phytoplankton_chlorophyll_ratio :: FT

    surface_intensity :: SI
    P_field :: P
    clock :: C
end

const TwoBandPhotosyntheticallyActiveRatiationField = Field{<:Center, <:Center, <:Center, <:TwoBandPhotosyntheticallyActiveRatiationOperand}

function TwoBandPhotosyntheticallyActiveRatiation(;grid, P, clock,
                                                   water_red_attenuation = 0.225, # 1/m
                                                   water_blue_attenuation = 0.0232, # 1/m
                                                   chlorophyll_red_attenuation = 0.037, # 1/(m * (mgChl/m³) ^ eʳ)
                                                   chlorophyll_blue_attenuation = 0.074, # 1/(m * (mgChl/m³) ^ eᵇ)
                                                   chlorophyll_red_exponent = 0.629,
                                                   chlorophyll_blue_exponent = 0.674,
                                                   pigment_ratio = 0.7,
                                                   phytoplankton_chlorophyll_ratio = 1.31, # mgChl/mol N
                                                   surface_intensity = t -> 100.0)
    return CenterField(grid, operand=TwoBandPhotosyntheticallyActiveRatiationOperand(water_red_attenuation,
                                                                                     water_blue_attenuation,
                                                                                     chlorophyll_red_attenuation,
                                                                                     chlorophyll_blue_attenuation,
                                                                                     chlorophyll_red_exponent,
                                                                                     chlorophyll_blue_exponent,
                                                                                     pigment_ratio,
                                                                                     phytoplankton_chlorophyll_ratio,
                                                                                     surface_intensity,
                                                                                     P,
                                                                                     clock))
end

function compute!(par::TwoBandPhotosyntheticallyActiveRatiationField)
    arch = architecture(par.grid)
    event = launch!(arch, par.grid, :xy, update_TwoBandPhotosyntheticallyActiveRatiation!, par, par.grid, par.operand.clock.time)
    wait(event)
end

"""
Sugar kelp model of [Broch2012](@cite) and updated by [Broch2013](@cite), [Fossberg2018](@cite), and [Broch2019](@cite).

Prognostic properties
===============
* Area: A (dm²)
* Nitrogen reserve: N (gN/gSW)
* Carbon reserve: C (gC/gSW)

Tracer dependencies
==============
* Nitrates: NO₃ (mmol N/m³)
* Photosynthetically available radiation: PAR (einstein/m²/day)

Optionally:
* Ammonia: NH₄ (mmol N/m³)
""" 
module SLatissimaModel
using Roots, KernelAbstractions
using OceanBioME.Particles: BiogeochemicalParticles, get_node
using Oceananigans.Units
using Oceananigans: CPU, Center
using Oceananigans.Architectures: arch_array, device, architecture
using Oceananigans: NonhydrostaticModel, HydrostaticFreeSurfaceModel
using Oceananigans.Biogeochemistry: required_biogeochemical_tracers, biogeochemical_auxiliary_fields

using KernelAbstractions.Extras.LoopInfo: @unroll 
using Oceananigans.Operators: volume
using Oceananigans.Grids: AbstractGrid
using Oceananigans.Fields: fractional_indices, _interpolate, datatuple

import Adapt: adapt_structure
import Oceananigans.Biogeochemistry: update_tendencies!
import Oceananigans.LagrangianParticleTracking: update_particle_properties!, _advect_particles!

@inline no_dynamics(args...) = nothing


"""
    SLatissima(; architecture :: AR = CPU()
                 growth_rate_adjustement :: FT = 4.5
                 photosynthetic_efficiency :: FT = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60)
                 minimum_carbon_reserve :: FT = 0.01
                 structural_carbon :: FT = 0.2
                 exudation :: FT = 0.5
                 erosion :: FT = 0.22
                 saturation_irradiance :: FT = 90 * day/ (10 ^ 6)
                 structural_dry_weight_per_area :: FT = 0.5
                 structural_dry_to_wet_weight :: FT = 0.0785
                 carbon_reserve_per_carbon :: FT = 2.1213
                 nitrogen_reserve_per_nitrogen :: FT = 2.72
                 minimum_nitrogen_reserve :: FT = 0.0126
                 maximum_nitrogen_reserve :: FT = 0.0216
                 growth_adjustement_2 :: FT = 0.039 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve))
                 growth_adjustement_1 :: FT = 0.18 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve)) - growth_adjustement_2
                 maximum_specific_growth_rate :: FT = 0.18
                 structural_nitrogen :: FT = 0.0146
                 photosynthesis_at_ref_temp_1 :: FT = 1.22e-3 * 24
                 photosynthesis_at_ref_temp_2 :: FT = 1.3e-3 * 24
                 photosynthesis_ref_temp_1 :: FT = 285.0
                 photosynthesis_ref_temp_2 :: FT = 288.0
                 photoperiod_1 :: FT = 0.85
                 photoperiod_2 :: FT = 0.3
                 respiration_at_ref_temp_1 :: FT = 2.785e-4 * 24
                 respiration_at_ref_temp_2 :: FT = 5.429e-4 * 24
                 respiration_ref_temp_1 :: FT = 285.0
                 respiration_ref_temp_2 :: FT = 290.0
                 photosynthesis_arrhenius_temp :: FT = (1 / photosynthesis_ref_temp_1 - 1 / photosynthesis_ref_temp_2) ^ -1 * log(photosynthesis_at_ref_temp_2 / photosynthesis_at_ref_temp_1)
                 photosynthesis_low_temp :: FT = 271.0
                 photosynthesis_high_temp :: FT = 296.0
                 photosynthesis_high_arrhenius_temp :: FT = 1414.87
                 photosynthesis_low_arrhenius_temp :: FT = 4547.89
                 respiration_arrhenius_temp :: FT = (1 / respiration_ref_temp_1 - 1 / respiration_ref_temp_2) ^ -1 * log(respiration_at_ref_temp_2 / respiration_at_ref_temp_1)
                 current_speed_for_0p65_uptake :: FT = 0.03
                 nitrate_half_saturation :: FT = 4.0
                 ammonia_half_saturation :: FT = 1.3
                 maximum_nitrate_uptake :: FT = 10 * structural_dry_weight_per_area * 24 * 14 / (10^6)
                 maximum_ammonia_uptake :: FT = 12 * structural_dry_weight_per_area * 24 * 14 / (10^6)
                 current_1 :: FT = 0.72
                 current_2 :: FT = 0.28
                 current_3 :: FT = 0.045
                 respiration_reference_A :: FT = 1.11e-4 * 24
                 respiration_reference_B :: FT = 5.57e-5 * 24
                 exudation_redfield_ratio :: FT = Inf

                 pescribed_velocity :: U = 0.1
                 pescribed_temperature :: T = nothing
                 pescribed_salinity :: S = nothing

                 #position
                 x :: P = arch_array(architecture, [0.0])
                 y :: P = arch_array(architecture, zeros(Float64, length(x)))
                 z :: P = arch_array(architecture, zeros(Float64, length(x)))

                 #properties
                 A :: P = arch_array(architecture, ones(Float64, length(x)) * 30)
                 N :: P = arch_array(architecture, ones(Float64, length(x)) * 0.01)
                 C :: P = arch_array(architecture, ones(Float64, length(x)) * 0.1)

                 #feedback
                 nitrate_uptake :: P = arch_array(architecture, zeros(Float64, length(x)))
                 ammonia_uptake :: P = arch_array(architecture, zeros(Float64, length(x)))
                 primary_production :: P = arch_array(architecture, zeros(Float64, length(x)))
                 frond_exudation :: P = arch_array(architecture, zeros(Float64, length(x)))
                 nitrogen_erosion :: P = arch_array(architecture, zeros(Float64, length(x)))
                 carbon_erosion :: P = arch_array(architecture, zeros(Float64, length(x)))

                 custom_dynamics :: F = no_dynamics

                 scalefactor :: FT = 1.0
                 latitude :: FT = 57.5)

Keywork Arguments
===================

- `architecture`: the architecture to adapt arrays to
- `growth_rate_adjustement`, ..., `exudation_redfield_ratio`: parameter values
- `pescribed_velocity`, `pescribed_temperature` and `pescribed_salinity`: functions for the relative velocity, temperature and salinity in the form `f(x, y, z, t)`
- `x`,`y` and `z`: positions of the particles
- `A`, `N`, and `C`: area, nitrogen, and carbon reserves
- `nitrate_uptake` ... `carbon_erosion`: diagnostic values coupled to tracer fields
- `custom_dynamics`: place to add any function of form `f!(particles, model, bgc, Δt)`
- `scalefactor`: scalar scaling for tracer coupling
- `latitude`: model latitude for seasonal growth modulation
"""
Base.@kwdef struct SLatissima{AR, FT, U, T, S, P, F} <: BiogeochemicalParticles
    architecture :: AR = CPU()
    growth_rate_adjustement :: FT = 4.5
    photosynthetic_efficiency :: FT = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60)
    minimum_carbon_reserve :: FT = 0.01
    structural_carbon :: FT = 0.2
    exudation :: FT = 0.5
    erosion :: FT = 0.22
    saturation_irradiance :: FT = 90 * day/ (10 ^ 6)
    structural_dry_weight_per_area :: FT = 0.5
    structural_dry_to_wet_weight :: FT = 0.0785
    carbon_reserve_per_carbon :: FT = 2.1213
    nitrogen_reserve_per_nitrogen :: FT = 2.72
    minimum_nitrogen_reserve :: FT = 0.0126
    maximum_nitrogen_reserve :: FT = 0.0216
    growth_adjustement_2 :: FT = 0.039 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve))
    growth_adjustement_1 :: FT = 0.18 / (2 * (1 - minimum_nitrogen_reserve / maximum_nitrogen_reserve)) - growth_adjustement_2
    maximum_specific_growth_rate :: FT = 0.18
    structural_nitrogen :: FT = 0.0146
    photosynthesis_at_ref_temp_1 :: FT = 1.22e-3 * 24
    photosynthesis_at_ref_temp_2 :: FT = 1.3e-3 * 24
    photosynthesis_ref_temp_1 :: FT = 285.0
    photosynthesis_ref_temp_2 :: FT = 288.0
    photoperiod_1 :: FT = 0.85
    photoperiod_2 :: FT = 0.3
    respiration_at_ref_temp_1 :: FT = 2.785e-4 * 24
    respiration_at_ref_temp_2 :: FT = 5.429e-4 * 24
    respiration_ref_temp_1 :: FT = 285.0
    respiration_ref_temp_2 :: FT = 290.0
    photosynthesis_arrhenius_temp :: FT = (1 / photosynthesis_ref_temp_1 - 1 / photosynthesis_ref_temp_2) ^ -1 * log(photosynthesis_at_ref_temp_2 / photosynthesis_at_ref_temp_1)
    photosynthesis_low_temp :: FT = 271.0
    photosynthesis_high_temp :: FT = 296.0
    photosynthesis_high_arrhenius_temp :: FT = 1414.87
    photosynthesis_low_arrhenius_temp :: FT = 4547.89
    respiration_arrhenius_temp :: FT = (1 / respiration_ref_temp_1 - 1 / respiration_ref_temp_2) ^ -1 * log(respiration_at_ref_temp_2 / respiration_at_ref_temp_1)
    current_speed_for_0p65_uptake :: FT = 0.03
    nitrate_half_saturation :: FT = 4.0
    ammonia_half_saturation :: FT = 1.3
    maximum_nitrate_uptake :: FT = 10 * structural_dry_weight_per_area * 24 * 14 / (10^6)
    maximum_ammonia_uptake :: FT = 12 * structural_dry_weight_per_area * 24 * 14 / (10^6)
    current_1 :: FT = 0.72
    current_2 :: FT = 0.28
    current_3 :: FT = 0.045
    respiration_reference_A :: FT = 1.11e-4 * 24
    respiration_reference_B :: FT = 5.57e-5 * 24
    exudation_redfield_ratio :: FT = Inf

    pescribed_velocity :: U = 0.1
    pescribed_temperature :: T = nothing
    pescribed_salinity :: S = nothing

    #position
    x :: P = arch_array(architecture, [0.0])
    y :: P = arch_array(architecture, zeros(Float64, length(x)))
    z :: P = arch_array(architecture, zeros(Float64, length(x)))

    #properties
    A :: P = arch_array(architecture, ones(Float64, length(x)) * 30)
    N :: P = arch_array(architecture, ones(Float64, length(x)) * 0.01)
    C :: P = arch_array(architecture, ones(Float64, length(x)) * 0.1)

    #feedback
    nitrate_uptake :: P = arch_array(architecture, zeros(Float64, length(x)))
    ammonia_uptake :: P = arch_array(architecture, zeros(Float64, length(x)))
    primary_production :: P = arch_array(architecture, zeros(Float64, length(x)))
    frond_exudation :: P = arch_array(architecture, zeros(Float64, length(x)))
    nitrogen_erosion :: P = arch_array(architecture, zeros(Float64, length(x)))
    carbon_erosion :: P = arch_array(architecture, zeros(Float64, length(x)))

    custom_dynamics :: F = no_dynamics

    scalefactor :: FT = 1.0
    latitude :: FT = 57.5
end

# for some reason this doesn't get called, 
# but if you manually call `adapt_structure(CuArray, particles::SLatissima)` 
# it returns a GPU friendly version. To automagically overcome this I'm `arch_array`ing 
# above but that does necessitate passing the grid to the particles
# weird thing is this is definitly getting called a some point, but not with the correct `to`.
adapt_structure(to, kelp::SLatissima) = SLatissima(kelp.architecture,
                                                   kelp.growth_rate_adjustement, 
                                                   kelp.photosynthetic_efficiency,
                                                   kelp.minimum_carbon_reserve,
                                                   kelp.structural_carbon,
                                                   kelp.exudation,
                                                   kelp.erosion,
                                                   kelp.saturation_irradiance,
                                                   kelp.structural_dry_weight_per_area,
                                                   kelp.structural_dry_to_wet_weight,
                                                   kelp.carbon_reserve_per_carbon,
                                                   kelp.nitrogen_reserve_per_nitrogen,
                                                   kelp.minimum_nitrogen_reserve,
                                                   kelp.maximum_nitrogen_reserve,
                                                   kelp.growth_adjustement_2,
                                                   kelp.growth_adjustement_1,
                                                   kelp.maximum_specific_growth_rate,
                                                   kelp.structural_nitrogen,
                                                   kelp.photosynthesis_at_ref_temp_1,
                                                   kelp.photosynthesis_at_ref_temp_2,
                                                   kelp.photosynthesis_ref_temp_1,
                                                   kelp.photosynthesis_ref_temp_2,
                                                   kelp.photoperiod_1,
                                                   kelp.photoperiod_2,
                                                   kelp.respiration_at_ref_temp_1,
                                                   kelp.respiration_at_ref_temp_2,
                                                   kelp.respiration_ref_temp_1,
                                                   kelp.respiration_ref_temp_2,
                                                   kelp.photosynthesis_arrhenius_temp,
                                                   kelp.photosynthesis_low_temp,
                                                   kelp.photosynthesis_high_temp,
                                                   kelp.photosynthesis_high_arrhenius_temp,
                                                   kelp.photosynthesis_low_arrhenius_temp,
                                                   kelp.respiration_arrhenius_temp,
                                                   kelp.current_speed_for_0p65_uptake,
                                                   kelp.nitrate_half_saturation,
                                                   kelp.ammonia_half_saturation,
                                                   kelp.maximum_nitrate_uptake,
                                                   kelp.maximum_ammonia_uptake,
                                                   kelp.current_1,
                                                   kelp.current_2,
                                                   kelp.current_3,
                                                   kelp.respiration_reference_A,
                                                   kelp.respiration_reference_B,
                                                   kelp.exudation_redfield_ratio,
                                                   kelp.pescribed_velocity,
                                                   kelp.pescribed_temperature,
                                                   kelp.pescribed_salinity,
                                                   adapt_structure(to, kelp.x),
                                                   adapt_structure(to, kelp.y),
                                                   adapt_structure(to, kelp.z),
                                                   adapt_structure(to, kelp.A),
                                                   adapt_structure(to, kelp.N),
                                                   adapt_structure(to, kelp.C),
                                                   adapt_structure(to, kelp.nitrate_uptake),
                                                   adapt_structure(to, kelp.ammonia_uptake),
                                                   adapt_structure(to, kelp.primary_production),
                                                   adapt_structure(to, kelp.frond_exudation),
                                                   adapt_structure(to, kelp.nitrogen_erosion),
                                                   adapt_structure(to, kelp.carbon_erosion),
                                                   kelp.custom_dynamics,
                                                   kelp.scalefactor,
                                                   kelp.latitude)

function update_tendencies!(bgc, particles::SLatissima, model)
    num_particles = length(particles)
    workgroup = min(num_particles, 256)
    worksize = num_particles

    update_tracer_tendencies_kernal! = update_tracer_tendencies!(device(model.architecture), workgroup, worksize)
    update_tracer_tendencies_event = update_tracer_tendencies_kernal!(bgc, particles, model.timestepper.Gⁿ, model.grid)
    wait(update_tracer_tendencies_event)
end

@kernel function update_tracer_tendencies!(bgc, p, tendencies, grid::AbstractGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    idx = @index(Global)

    x = p.x[idx]
    y = p.y[idx]
    z = p.z[idx]

    bgc_tracers = required_biogeochemical_tracers(bgc)

    #nodes, weights, normfactor = @inbounds get_nearest_nodes(x, y, z, grid, (Center(), Center(), Center()))

    #@unroll for (n, weight) in enumerate(weights)
        # Reflect back on Bounded boundaries or wrap around for Periodic boundaries
    i, j, k = fractional_indices(x, y, z, (Center(), Center(), Center()), grid)

    # Convert fractional indices to unit cell coordinates 0 <= (ξ, η, ζ) <=1
    # and integer indices (with 0-based indexing).
    ξ, i = modf(i)
    η, j = modf(j)
    ζ, k = modf(k)

    i, j, k = (get_node(TX(), Int(ifelse(ξ < 0.5, i + 1, i + 2)), grid.Nx), get_node(TY(), Int(ifelse(η < 0.5, j + 1, j + 2)), grid.Ny), get_node(TZ(), Int(ifelse(ζ < 0.5, k + 1, k + 2)), grid.Nz))

    node_volume = volume(i, j, k, grid, Center(), Center(), Center())

    node_scalefactor = p.scalefactor / node_volume #* normfactor / (weight * node_volume)

    @inbounds begin
        tendencies.NO₃[i, j, k] -= node_scalefactor * p.nitrate_uptake[idx]
        
        if :NH₄ in bgc_tracers
            tendencies.NH₄[i, j, k] -= node_scalefactor * p.ammonia_uptake[idx]
        end

        if :DIC in bgc_tracers
            tendencies.DIC[i, j, k] -= node_scalefactor * p.primary_production[idx]
        end

        if :O₂ in bgc_tracers
            tendencies.O₂[i, j, k] += node_scalefactor * p.primary_production[idx]
        end

        if :DOM in bgc_tracers
            tendencies.DOM[i, j, k] += node_scalefactor * p.frond_exudation[idx] / p.exudation_redfield_ratio
        elseif :DON in bgc_tracers
            tendencies.DON[i, j, k] += node_scalefactor * p.frond_exudation[idx] / p.exudation_redfield_ratio
            tendencies.DOC[i, j, k] += node_scalefactor * p.frond_exudation[idx]
        end

        if :bPOM in bgc_tracers
            tendencies.bPOM[i, j, k] += node_scalefactor * p.nitrogen_erosion[idx]
        elseif :bPON in bgc_tracers
            tendencies.bPON[i, j, k] += node_scalefactor * p.nitrogen_erosion[idx]
            tendencies.bPOC[i, j, k] += node_scalefactor * p.carbon_erosion[idx]
        end
    end
end

function update_particle_properties!(particles::SLatissima, model, bgc, Δt)
    workgroup = min(length(particles), 256)
    worksize = length(particles)

    arch = architecture(model)

    # Advect particles
    advect_particles_kernel! = _advect_particles!(device(arch), workgroup, worksize)

    advect_particles_event = advect_particles_kernel!((x = particles.x, y = particles.y, z = particles.z), 
                                                      1.0, model.grid, Δt,
                                                      datatuple(model.velocities),
                                                      dependencies=Event(device(arch)))

    wait(device(arch), advect_particles_event)

    update_particle_properties_kernel! = _update_particle_properties!(device(arch), workgroup, worksize)

    update_particle_properties_event = update_particle_properties_kernel!(particles, bgc, model, Δt,
                                                                          dependencies=Event(device(arch)))

    wait(device(arch), update_particle_properties_event)

    particles.custom_dynamics(particles, model, bgc, Δt)
end

@kernel function _update_particle_properties!(p, bgc, model, Δt)
    idx = @index(Global)

    @inbounds begin
        x = p.x[idx]
        y = p.y[idx]
        z = p.z[idx]

        A = p.A[idx]
        N = p.N[idx]
        C = p.C[idx]
    end

    t = model.clock.time

    @inbounds if p.A[idx] > 0.0
        NO₃, NH₄, PAR, u, T, S = get_arguments(x, y, z, t, p, bgc, model)

        photo = photosynthesis(T, PAR, p)
        e = exudation(C, p)
        ν = erosion(A, p)

        fᶜ = f_curr(u, p)

        j_NO₃ = max(0.0, p.maximum_ammonia_uptake * fᶜ * ((p.maximum_nitrogen_reserve - N) /
                         (p.maximum_nitrogen_reserve - p.minimum_nitrogen_reserve)) * NO₃ / 
                         (p.nitrate_half_saturation + NO₃))

        j̃_NH₄ = max(0.0, p.maximum_ammonia_uptake * fᶜ * NH₄ / (p.ammonia_half_saturation + NH₄))

        μ_NH₄ = j̃_NH₄ / (p.structural_dry_weight_per_area * (N + p.structural_nitrogen))
        μ_NO₃ = 1 - p.minimum_nitrogen_reserve / N
        μ_C = 1 - p.minimum_carbon_reserve / C

        n = floor(Int, mod(t, 364days) / day)
        λ = normed_day_length_change(n, p.latitude)

        μ = @inbounds f_area(A, p) * 
                      f_temp(T, p) * 
                      f_photo(λ, p) * 
                      min(μ_C, max(μ_NO₃, μ_NH₄))

        j_NH₄ = min(j̃_NH₄, μ * p.structural_dry_weight_per_area * (N + p.structural_nitrogen))

        r = respiration(T, μ, j_NO₃ + j_NH₄, p)

        dA = (μ - ν) * A / day

        dN = ((j_NO₃ + j_NH₄ - photo * e * 14 / (12 * p.exudation_redfield_ratio)) / p.structural_dry_weight_per_area - 
               μ * (N + p.structural_nitrogen)) / day

        dC = ((photo * (1 - e) - r) / p.structural_dry_weight_per_area - 
              μ * (C + p.structural_carbon)) / day

        A_new = A + dA * Δt 
        N_new = N + dN * Δt 
        C_new = C + dC * Δt

        if C_new < p.minimum_carbon_reserve
            A_new *= (1 - (p.minimum_carbon_reserve - C) / p.structural_carbon)
            C_new = p.minimum_carbon_reserve
            N_new += p.structural_nitrogen * (p.minimum_carbon_reserve - C) / p.structural_carbon
        end
        
        if N_new < p.minimum_nitrogen_reserve
            A_new *= (1 - (p.minimum_nitrogen_reserve - N) / p.structural_nitrogen)
            N_new = p.minimum_nitrogen_reserve
            C_new += p.structural_carbon * (p.minimum_nitrogen_reserve - N) / p.structural_nitrogen
        end

        @inbounds begin
            p.primary_production[idx] = (photo - r) * A / (day * 12 * 0.001) #gC/dm^2/hr to mmol C/s

            p.frond_exudation[idx] = e * photo * A / (day * 12 * 0.001)#mmol C/s

            p.nitrogen_erosion[idx] = ν * p.structural_dry_weight_per_area * A * (p.structural_nitrogen + N) / (day * 14 * 0.001)#1/day to mmol N/s

            p.carbon_erosion[idx] = ν * p.structural_dry_weight_per_area * A * (p.structural_carbon + C) / (day * 12 * 0.001)#1/hr to mmol C/s

            p.nitrate_uptake[idx] = j_NO₃ * A / (day * 14 * 0.001)#gN/dm^2/hr to mmol N/s

            p.ammonia_uptake[idx] = j_NH₄ * A / (day * 14 * 0.001)#gN/dm^2/hr to mmol N/s

            p.A[idx] = A_new
            p.N[idx] = N_new
            p.C[idx] = C_new
        end
    end
end


@inline function photosynthesis(T, PAR, p)
    Tk = T + 273.15

    pₘₐₓ = p.photosynthesis_at_ref_temp_1 * 
            exp(p.photosynthesis_arrhenius_temp / p.photosynthesis_ref_temp_1 -
                p.photosynthesis_arrhenius_temp / Tk) / 
            (1 + exp(p.photosynthesis_low_arrhenius_temp / Tk -
                     p.photosynthesis_low_arrhenius_temp / p.photosynthesis_low_temp) 
               + exp(p.photosynthesis_high_arrhenius_temp / p.photosynthesis_high_temp -
                     p.photosynthesis_high_arrhenius_temp / Tk))

    β = find_zero(β_func, (0, 0.1), Bisection(); p = (; pₘₐₓ, α = p.photosynthetic_efficiency, I_sat = p.saturation_irradiance))

    pₛ = p.photosynthetic_efficiency * p.saturation_irradiance / log(1 + p.photosynthetic_efficiency / β)

    return pₛ * (1 - exp(- p.photosynthetic_efficiency * PAR / pₛ)) * exp(-β * PAR / pₛ) 
end

@inline β_func(x, p) = p.pₘₐₓ - (p.α * p.I_sat / log(1 + p.α / x)) * 
                                (p.α / (p.α + x)) * (x / (p.α + x))^(x / p.α)

@inline exudation(C, p) = 1 - exp(p.exudation * (p.minimum_carbon_reserve - C))

@inline erosion(A, p) = 1e-6 * exp(p.erosion * A) / (1 + 1e-6 * (exp(p.erosion * A) - 1))

@inline respiration(T, μ, j, p) = (p.respiration_reference_A * (μ / p.maximum_specific_growth_rate + 
                                                                j / (p.maximum_nitrate_uptake + p.maximum_ammonia_uptake)) +
                                   p.respiration_reference_B) * 
                                  exp(p.respiration_arrhenius_temp / p.respiration_ref_temp_1 - p.respiration_arrhenius_temp / (T + 273.15))

#####
##### Growth limitation
#####

@inline f_curr(u, p) = p.current_1 * (1 - exp(-u / p.current_3)) + p.current_2

@inline f_area(a, p) = p.growth_adjustement_1 * exp(-(a / p.growth_rate_adjustement)^2) + p.growth_adjustement_2

@inline function f_temp(T, p)
    # should probably parameterise these limits
    if -1.8 <= T < 10 
        return 0.08 * T + 0.2
    elseif 10 <= T <= 15
        return 1
    elseif 15 < T <= 19
        return 19 / 4 - T / 4
    elseif T > 19
        return 0
    else
        return 0
    end
end

@inline f_photo(λ, p) = p.photoperiod_1 * (1 + sign(λ) * abs(λ) ^ 0.5) + p.photoperiod_2

@inline function day_length(n, ϕ)
    n -= 171
    M = mod((356.5291 + 0.98560028 * n), 360)
    C = 1.9148 * sin(M * π / 180) + 0.02 * sin(2 * M * π / 180) + 0.0003 * sin(3 * M * π / 180)
    λ = mod(M + C + 180 + 102.9372, 360)
    δ = asin(sin(λ * π / 180) * sin(23.44 * π / 180))
    ω = (sin(-0.83 * π / 180) * sin(ϕ * π / 180) * sin(δ))/(cos(ϕ * π / 180) * cos(δ))
    return ω / 180
end

@inline normed_day_length_change(n, ϕ) = (day_length(n, ϕ) - day_length(n - 1, ϕ)) / (day_length(76, ϕ) - day_length(75, ϕ))


@inline function get_arguments(x, y, z, t, particles, bgc, model)
    bgc_tracers = required_biogeochemical_tracers(bgc)
    auxiliary_fields = biogeochemical_auxiliary_fields(bgc)

    i, j, k = fractional_indices(x, y, z, (Center(), Center(), Center()), model.grid)

    ξ, i = modf(i)
    η, j = modf(j)
    ζ, k = modf(k)

    # TODO: ADD ALIASING/RENAMING OF TRACERS SO WE CAN USE E.G. N IN STEAD OF NO3

    NO₃ = _interpolate(model.tracers.NO₃, ξ, η, ζ, Int(i+1), Int(j+1), Int(k+1))
    PAR = _interpolate(auxiliary_fields.PAR, ξ, η, ζ, Int(i+1), Int(j+1), Int(k+1)) * day / (3.99e-10 * 545e12) # W / m² / s to einstein / m² / day

    if :NH₄ in bgc_tracers
        NH₄ = _interpolate(model.tracers.NH₄, ξ, η, ζ, Int(i+1), Int(j+1), Int(k+1))
    else
        NH₄ = 0.0
    end

    if isnothing(particles.pescribed_velocity)
        u =  sqrt(_interpolate(model.velocities.u, ξ, η, ζ, Int(i+1), Int(j+1), Int(k+1)) .^ 2 + 
                  _interpolate(model.velocities.v, ξ, η, ζ, Int(i+1), Int(j+1), Int(k+1)) .^ 2 + 
                  _interpolate(model.velocities.w, ξ, η, ζ, Int(i+1), Int(j+1), Int(k+1)) .^ 2)
    elseif isa(particles.pescribed_velocity, Number)
        u = particles.pescribed_velocity
    else
        u = particles.pescribed_velocity(x, y, z, t)
    end

    if isnothing(particles.pescribed_temperature)
        T = _interpolate(model.tracers.T, ξ, η, ζ, Int(i+1), Int(j+1), Int(k+1))
    elseif isa(particles.pescribed_temperature, Number)
        T = particles.pescribed_temperature
    else
        T = particles.pescribed_temperature(x, y, z, t)
    end

    #=if isnothing(particles.pescribed_salinity)
        S = _interpolate(model.tracers.S, ξ, η, ζ, Int(i+1), Int(j+1), Int(k+1))
    else
        S = particles.pescribed_salinity(x, y, z, t)
    end=#
    S = 35.0 # we currently don't use S for anything

    return NO₃, NH₄, PAR, u, T, S
end
end

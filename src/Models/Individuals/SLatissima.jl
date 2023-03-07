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

Tracer coupling
===========
* j_NO₃ and j_NH₄: uptake of nitrates and Ammonia
* j_DIC and j_OXY: update of disolved inorganic carbon and release/uptake of oxygen
* e: exudation of disolved organic matter
* νⁿ and νᶜ: loss of mass as large particles (nitrogen and carbon content)

Notes
=====
SW is the structural weight and is equal to K_A*A
"""
module SLatissima
using StructArrays, Roots
using OceanBioME.Particles: ActiveLagrangianParticles
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

import Adapt: adapt_structure

#####
##### Photosynthesis
#####

@inline β_func(x, params) = params.p_max - (params.α * params.I_sat / log(1 + params.α / x)) * (params.α / (params.α + x)) * (x / (params.α + x))^(x / params.α)

@inline function _p(temp, irr, params)
    p_max = params.P_1 * exp(params.T_AP / params.T_P1 - params.T_AP / (temp + 273.15)) / (1 + exp(params.T_APL / (temp + 273.15) - params.T_APL / params.T_PL) + exp(params.T_APH / params.T_PH - params.T_APH / (temp + 273.15)))
    β = find_zero(β_func, (0, 0.1), Bisection(); p=merge(params, (; p_max)))
    p_s = params.α * params.I_sat / log(1 + params.α / β)
    return p_s * (1 - exp(-params.α * irr / p_s)) * exp(-β * irr / p_s) 
end

#####
##### Mass loss
#####

@inline _e(c, params) = 1 - exp(params.γ * (params.C_min - c))

@inline _ν(a, params) = 1e-6 * exp(params.ϵ * a) / (1 + 1e-6 * (exp(params.ϵ * a) - 1))

#####
##### Growth limitation
#####

@inline f_curr(u, params) = params.uₐ*(1-exp(-u/params.u₀))+params.uᵦ

@inline f_area(a, params) = params.m_1 * exp(-(a / params.A_0)^2) + params.m_2

@inline function f_temp(temp, params)
    if -1.8 <= temp < 10 
        return 0.08 * temp + 0.2
    elseif 10 <= temp <= 15
        return 1
    elseif 15 < temp <= 19
        return 19 / 4 - temp / 4
    elseif temp > 19
        return 0
    else
        return 0
    end
end

@inline f_photo(λ, params) = params.a_1 * (1 + sign(λ) * abs(λ)^.5) + params.a_2

#####
##### Respiration
#####
@inline _r(temp, μ, j, params) = (params.R_A * (μ / params.μ_max + j / params.J_max) + params.R_B) * exp(params.T_AR / params.T_R1 - params.T_AR / (temp + 273.15))

#####
##### Growth equations
#####

function equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, params, Δt::AbstractFloat)
    if !iszero(A)
        irr /= 3.99e-10*545e12/(1day) #W/m²/s to einstein/m²/day
        p = _p(T, irr, params)
        e = _e(C, params)
        ν  = _ν(A, params)

        j_NO₃ = max(0.0, params.j_NO₃_max * f_curr(u, params) * ((params.N_max - N) / (params.N_max - params.N_min)) * NO₃ / (params.k_NO₃ + NO₃))
        j̃_NH₄ = max(0.0, params.j_NH₄_max * f_curr(u, params) * NH₄ / (params.k_NH₄ + NH₄))

        μ_NH₄ = j̃_NH₄ / (params.K_A * (N + params.N_struct))
        μ_NO₃ = 1 - params.N_min / N
        μ_C = 1 - params.C_min / C

        μ = @inbounds f_area(A, params) * f_temp(T, params) * f_photo(params.λ[1 + floor(Int, mod(t, 364days) / day)], params) * min(μ_C, max(μ_NO₃, μ_NH₄))

        j_NH₄ = min(j̃_NH₄, μ * params.K_A * (N + params.N_struct))

        r = _r(T, μ, j_NO₃ + j_NH₄, params)

        dA = (μ - ν) * A / (60*60*24)
        dN = ((j_NO₃ + j_NH₄ - p * e * 14/(12 * params.exudation_redfield_ratio)) / params.K_A - μ * (N + params.N_struct)) / (60 * 60 * 24)
        dC = ((p * (1 - e) - r) / params.K_A - μ * (C + params.C_struct)) / (60 * 60 * 24)

        A_new = A+dA*Δt 
        N_new = N+dN*Δt 
        C_new = C+dC*Δt

        if C_new < params.C_min
            A_new *= (1 - (params.C_min - C) / params.C_struct)
            C_new = params.C_min
            N_new += params.N_struct * (params.C_min - C) / params.C_struct
        end
        
        if N_new < params.N_min
            A_new *= (1 - (params.N_min - N) / params.N_struct)
            N_new = params.N_min
            C_new += params.C_struct * (params.N_min - N) / params.N_struct
        end

        pp = (p - r) * A / (60 * 60 * 24 * 12 * 0.001) #gC/dm^2/hr to mmol C/s
        e *= p * A / (60 * 60 * 24 * 12 * 0.001)#mmol C/s
        νⁿ = ν * params.K_A * A * (params.N_struct + N) / (60 * 60 * 24 * 14 * 0.001)#1/hr to mmol N/s
        νᶜ = ν * params.K_A * A * (params.C_struct + C) / (60 * 60 * 24 * 12 * 0.001)#1/hr to mmol C/s
        j_NO₃ *= A / (60 * 60 * 24 * 14 * 0.001)#gN/dm^2/hr to mmol N/s
        j_NH₄ *= A / (60 * 60 * 24 * 14 * 0.001)#gN/dm^2/hr to mmol N/s

        du, dv, dw = 0.0, 0.0, 0.0
    else
        A_new, N_new, C_new, j_NO₃, j_NH₄, pp, e, νⁿ, νᶜ, du, dv, dw = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end
    return (A = A_new, N = N_new, C = C_new, j_NO₃ = -j_NO₃, j_NH₄ = -j_NH₄, j_DIC = -pp, j_OXY=pp, eᶜ = e, eⁿ = e / params.exudation_redfield_ratio, νⁿ = νⁿ, νᶜ = νᶜ, u = du, v = dv, w = dw)
end

#fixed urel, T and S functions
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, irr::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, params.urel, params, Δt)
#fixed urel, T and S fields
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, params.urel, params, Δt)
#tracked u, T and S functions can not be done (like this) bcs same number of variables as main function
 #equations(x::AbstractFloat, y::AbstractFloat, z::Abstractparams.T(x, y, z, t)Float, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, v::AbstractFloat, w::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, sqrt(u^2+v^2+w^2), params, Δt)
#tracked u, T and S fields
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, v::AbstractFloat, w::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, T, S, irr, sqrt(u^2+v^2+w^2), params, Δt)

#####
##### Parameters (some diagnostic so have to define before)
#####

N_min = 0.0126
N_max = 0.0216
m_2 = 0.039/(2*(1-N_min/N_max))
T_P1 = 285
T_P2 = 288
P_1 = 1.22e-3 * 24
P_2 = 1.3e-3 * 24
T_R1 = 285
T_R2 = 290
R_1 =  2.785e-4 * 24
R_2 = 5.429e-4 * 24
K_A = 0.5

const defaults = (
    A_0 = 4.5,#6,# Growth rate adjustment parameter
    α = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60),#3.75e-5 * 24 * 10^6 / (24 * 60 * 60),# photosynthetic efficiency 
    C_min = 0.01,# Minimal carbon reserve
    C_struct = 0.2,# Amount of carbon per unit dry weight of structural mass
    γ = 0.5,# Exudation parameter
    ϵ = 0.22,# Frond erosion parameter
    I_sat = 90 * 24*60*60/(10^6),#200 * 24 * 60 * 60 / (10^6),# Irradiance for maximal photosynthesis
    J_max = 1.4e-4 * 24,# Maximal nitrate uptake (gN/dm^2/h converted to gN/dm^2/day)
    K_A = K_A,#0.6,# Structural dry weight per unit area
    K_DW = 0.0785,# Dry weight to wet weight ratio of structural mass
    K_C = 2.1213,# Mass of carbon reserves per gram carbon
    K_N = 2.72,# Mass of nitrogen reserves per gram nitrogen
    N_min = N_min,#0.01,# Minimal nitrogen reserve
    N_max = N_max,#0.022,# Maximal nitrogen reserve
    m_2 = m_2,#0.036,# 0.03,# Growth rate adjustment parameter
    m_1 = 0.18/(2*(1-N_min/N_max))-m_2,# 0.1085,# Growth rate adjustment parameter
    μ_max = 0.18,# Maximal area specific growth ratio
    N_struct = 0.0146,#0.01,# Amount of nitrogen per unit dry weight of structural mass
    P_1 = P_1,# Maximal photosynthetic rate at T = T?P1K converted to day^-1
    P_2 = P_2,# Maximal photosynthetic rate at T = T?P2K converted to day^-1
    a_1 = 0.85,# Photoperiod parameter
    a_2 = 0.3,# Photoperiod parameter
    R_1 =  R_1,# Respiration rate at T = TR1,
    R_2 = R_2,# Respiration rate at T = TR2,
    T_R1 = T_R1,# Reference temperature for respiration (K)
    T_R2 = T_R2,# Reference temperature for respiration (K)
    T_P1 = T_P1,# Reference temperature for photosynthesis (K)
    T_P2 = T_P2,# Reference temperature for photosynthesis (K)
    T_AP = (1 / T_P1 - 1 / T_P2)^(-1) * log(P_2 / P_1),# 1694.4,# Arrhenius temperature for photosynthesis (K)
    T_PL = 271,
    T_PH = 296,
    T_APH = 1414.87,# 1516.7,#25924,# Arrhenius temperature for photosynthesis at high end of range (K)
    T_APL = 4547.89,# 4388.5,#27774,# Arrhenius temperature for photosynthesis at low end of range (K)
    T_AR = (1 / T_R1 - 1 / T_R2)^(-1) * log(R_2 / R_1),# Arrhenius temperature for respiration (K)
    U_0p65 = 0.03,# Current speed at which J = 0.65Jmax (m/s)
    k_NO₃ = 4,# Nitrate uptake half saturation constant, converted to mmol/m^3
    k_NH₄ = 1.3, #Fossberg 2018

    j_NO₃_max = 10 * K_A * 24 * 14 / (10^6), #Ahn 1998
    j_NH₄_max = 12 * K_A * 24 * 14 / (10^6), #Ahn 1998

    uₐ = 0.72, #broch 2019
    uᵦ = 0.28,
    u₀ = 0.045,

    R_A=1.11e-4*24,
    R_B=5.57e-5*24,

    exudation_redfield_ratio = Inf,
)

#####
##### Seasonality function
#####

function generate_seasonality(lat)
    θ = 0.2163108 .+ 2 .* atan.(0.9671396 .* tan.(0.00860 .* ([1:365;] .- 186)))
    δ = asin.(0.39795 .* cos.(θ))
    p = 0.8333
    d = 24 .- (24 / pi) .* acos.((sin(p * pi / 180) .+ sin(lat * pi / 180) .* sin.(δ)) ./ (cos(lat * pi / 180) .* cos.(δ)))
    λ = diff(d)
    push!(λ, λ[end])
    return λ ./ findmax(λ)[1]
end

#####
##### Model definition for use in OceanBioME
#####

struct Particle{FT}
    #position
    x :: FT
    y :: FT
    z :: FT
    #velocity
    u :: FT
    v :: FT
    w :: FT
    #properties
    A :: FT
    N :: FT
    C :: FT
    #feedback
    j_NO₃ :: FT
    j_NH₄ :: FT
    j_DIC :: FT
    j_OXY :: FT
    eⁿ :: FT
    eᶜ :: FT
    νⁿ :: FT
    νᶜ :: FT
    #tracked fields
    NO₃  :: FT
    NH₄  :: FT
    PAR :: FT
    T :: FT
    S :: FT
end

default_tracers = (NO₃ = :NO₃, PAR=:PAR) # tracer = property
optional_tracer_dependencies = (NH₄ = :NH₄, )
default_coupling = (NO₃ = :j_NO₃, )    
optional_tracer_coupling = (NH₄ = :j_NH₄, DIC = :j_DIC, DON = :eⁿ, DOC = :eᶜ, bPON = :νⁿ, bPOC = :νᶜ, O₂ = :j_OXY)

# Expands initial conditions to fill particles
function defineparticles(initials, n)
    x̄₀ = []
    for var in [:x₀, :y₀, :z₀, :A₀, :N₀, :C₀]
        vals=getproperty(initials, var)
        if isa(vals, AbstractFloat)
            push!(x̄₀, repeat([vals], n))
        elseif (isa(vals, Vector) && length(vals) == n)
            push!(x̄₀, vals)
        else
            throw(ArgumentError("Invalid initial values given for $var, must be a single number or vector of length n"))
        end
    end
    return StructArray{Particle}((x̄₀[1], x̄₀[2], x̄₀[3], zeros(n), zeros(n), zeros(n), x̄₀[4], x̄₀[5], x̄₀[6], [zeros(n) for i in 1:13]...))
end

@inline no_dynamics(args...) = nothing

"""
    SLatissima.setup(n, x₀, y₀, z₀, A₀, N₀, C₀, latitude;
                     scalefactor = 1.0, 
                     T=nothing, 
                     S=nothing, 
                     urel=nothing, 
                     paramset=defaults, 
                     custom_dynamics=no_dynamics, 
                     optional_tracers=(),
                     tracer_names=NamedTuple())

Returns Oceananigans `LagrangianParticles` which will evolve and interact with biogeochemistry as sugar kelp

Keyword arguments
=============

    - `n`: number of particles
    - `x₀`, `y₀`, `z₀`: intial position, should be array of length `n`
    - `A₀`, `N₀`, `C₀`: Initial area/nitrogen reserve/carbon reserve, can be arrays of length `n` or single values
    - `latitude`: latitude of growth (used for seasonality parameterisation)
    - `scalefactor`: multiplier on particle uptake/release of tracers, can be though of as the number of kelp frond represented by each particle
    - `T` and `S`: functions of form `func(x, y, z, t)` for the temperature and salinity, if not defined then model will look for tracer fields (useful if physics are resolved)
    - `urel`: relative water velocity (single number), if this is not specified then actual local relative water velocity will be used (moderates nutrient uptake)
    - `paramset`: parameters for growth model
    - `custom_dynamics`: any other dynamics you wish to be run on the particles as they evolve, should be of the form `dynamics!(particles, model, Δt)`
    - `optional_tracers`: Tuple of optional tracers to couple with
    - `tracer_names`: rename coupled tracers with NamedTuple with keys as new name and values as the property this feeds/depends on (e.g. `(N = :NO₃, )`)
"""

function setup(; n, x₀, y₀, z₀, A₀, N₀, C₀, latitude,
                 scalefactor = 1.0, 
                 T = nothing, 
                 S = nothing, 
                 urel = nothing, 
                 paramset = defaults, 
                 custom_dynamics = no_dynamics, 
                 optional_tracers = (),
                 tracer_names = NamedTuple())

    # this is a mess, we're not even using salinity at the moment
    if (!isnothing(T) && !isnothing(S) && isnothing(urel))
        throw(ArgumentError("T and S functions with tracked velocity fields not currently implimented"))
    elseif ((isnothing(T) && !isnothing(S)) | (!isnothing(T) && isnothing(S)))
        throw(ArgumentError("T and S must both be functions or both be tracked fields"))
    end

    # fills out particles in case we weren't given arrays
    particles = defineparticles((x₀=x₀, y₀=y₀, z₀=z₀, A₀=A₀, N₀=N₀, C₀=C₀), n)

    property_dependencies = (:A, :N, :C, :NO₃, :NH₄, :PAR)
    λ_arr = generate_seasonality(latitude)
    parameters = merge(paramset, (λ=λ_arr, ))
    diagnostic_properties = (:A, :N, :C, :j_NO₃, :j_NH₄, :j_DIC, :j_OXY, :eⁿ, :eᶜ, :νⁿ, :νᶜ) # all diagnostic for the sake of enforcing C limit

    if isnothing(T) property_dependencies = (property_dependencies..., :T, :S) else parameters = merge(parameters, (T=T, S=S)) end
    if isnothing(urel) property_dependencies = (property_dependencies..., :u, :v, :w) else parameters = merge(parameters, (urel = urel, )) end

    prognostic_properties = (:u, :v, :w) #for when I impliment dynamics

    tracers = default_tracers
    coupled = default_coupling
    for tracer in optional_tracers
        if tracer in optional_tracer_dependencies
            tracers = merge(tracers, NamedTuple{(tracer, )}((getproperty(optional_tracer_dependencies, tracer), )))
        end
        if tracer in keys(optional_tracer_coupling)
            coupled = merge(coupled, NamedTuple{(tracer, )}((getproperty(optional_tracer_coupling, tracer), )))
        end
    end

    # over writing tracer names with chosen names, would be a lot easier if we did property = tracer
    tracers_inverse = NamedTuple{values(tracers)}(keys(tracers))
    coupled_inverse = NamedTuple{values(coupled)}(keys(coupled))

    for (new_name, tracer) in pairs(tracer_names)
        if tracer in keys(tracers)
            tracer_inverse[tracers[tracer]] = new_name
        end
        
        if tracer in keys(coupled)
            coupled_inverse[coupled[tracer]] = new_name
        end
    end

    tracers = NamedTuple{values(tracers_inverse)}(keys(tracers_inverse))
    coupled = NamedTuple{values(coupled_inverse)}(keys(coupled_inverse))


    return ActiveLagrangianParticles(particles;
                                     equation = equations, 
                                     equation_arguments = property_dependencies, 
                                     equation_parameters = parameters, 
                                     prognostic = prognostic_properties, 
                                     diagnostic = diagnostic_properties, 
                                     tracked_fields = tracers, 
                                     coupled_fields = coupled,
                                     scalefactor = scalefactor, 
                                     custom_dynamics = custom_dynamics)
end
end

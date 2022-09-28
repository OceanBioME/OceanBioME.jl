"
Model of Sacharina Latissima from Broch and Slagstad 2012, updated with parameters of Broch 2013.
TODO: Also extended to depend on salinity and current by Broch 2019.

Extending to consider uptake of NO₃ vs NH₄ by Fossberg 2018, with parameter values from Ahn 1998

References:
Ahn, O., Petrell, R. J., & Harrison, P. J. (1998). Ammonium and nitrate uptake by Laminaria saccharina and Nereocystis luetkeana originating from a salmon sea cage farm. In Journal of Applied Phycology (Vol. 10, pp. 333–340)
Broch, OJ. and Slagstad, D., 2012. Modelling seasonal growth and composition of the kelp Saccharina latissima. Journal of Applied Phycology, 24(4), pp.759-776.
Broch, OJ., Ellingsen, I. H., Forbord, S., Wang, X., Volent, Z., Alver, M. O., Handå, A., Andresen, K., Slagstad, D., Reitan, K. I., Olsen, Y., & Skjermo, J. (2013). Modelling the cultivation and bioremediation potential of the kelp Saccharina latissima in close proximity to an exposed salmon farm in Norway. Aquaculture Environment Interactions, 4(2), 187–206. https://doi.org/10.3354/aei00080 
Broch OJ., Alver MO, Bekkby T, Gundersen H, Forbord S, Handå A, Skjermo J and Hancke K (2019) The Kelp Cultivation Potential in Coastal and Offshore Regions of Norway. Front. Mar. Sci. 5:529. doi: 10.3389/fmars.2018.00529
Fossberg, J., Forbord, S., Broch, O. J., Malzahn, A. M., Jansen, H., Handå, A., Førde, H., Bergvik, M., Fleddum, A. L., Skjermo, J., & Olsen, Y. (2018). The potential for upscaling kelp (Saccharina latissima) cultivation in salmon-driven integrated multi-trophic aquaculture (IMTA). Frontiers in Marine Science, 9(NOV). https://doi.org/10.3389/fmars.2018.00418 
"
module SLatissima
using StructArrays, Roots
using OceanBioME: Particles
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

@inline f_curr(u, params) = params.uₐ*(1-exp(-u/params.u₀))+params.uᵦ

#Restated from SugarKelp.jl to reduce dependency maintinance complexity
@inline β_func(x, params) = params.p_max - (params.α * params.I_sat / log(1 + params.α / x)) * (params.α / (params.α + x)) * (x / (params.α + x))^(x / params.α)
@inline function _p(temp, irr, params)
    p_max = params.P_1 * exp(params.T_AP / params.T_P1 - params.T_AP / (temp + 273.15)) / (1 + exp(params.T_APL / (temp + 273.15) - params.T_APL / params.T_PL) + exp(params.T_APH / params.T_PH - params.T_APH / (temp + 273.15)))
    β = find_zero(β_func, (0, 0.1), Bisection(); p=merge(params, (; p_max)))
    p_s = params.α * params.I_sat / log(1 + params.α / β)
    return p_s * (1 - exp(-params.α * irr / p_s)) * exp(-β * irr / p_s) 
end

@inline _e(c, params) = 1 - exp(params.γ * (params.C_min - c))
@inline _ν(a, params) = 1e-6 * exp(params.ϵ * a) / (1 + 1e-6 * (exp(params.ϵ * a) - 1))
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

@inline _r(temp, μ, j, params) = (params.R_A * (μ / params.μ_max + j / params.J_max) + params.R_B) * exp(params.T_AR / params.T_R1 - params.T_AR / (temp + 273.15))

function equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, params, Δt::AbstractFloat)
    if !iszero(A)
        irr /= 3.99e-10*545e12/(1day) #W/m²/s to einstein/m²/day
        p = _p(T, irr, params)
        e = _e(C, params)
        ν  = _ν(A, params)

        j_NO₃ = params.j_NO₃_max*f_curr(u, params)*((params.N_max-N)/(params.N_max-params.N_min))*NO₃/(params.k_NO₃+NO₃)
        j̃_NH₄ = params.j_NH₄_max*f_curr(u, params)*NH₄/(params.k_NH₄+NH₄)
        μ_NH₄ = j̃_NH₄/(params.K_A*(N+params.N_struct))
        μ_NO₃ = 1 - params.N_min/N
        μ_C = 1 - params.C_min/C

        μ = @inbounds f_area(A, params)*f_temp(T, params)*f_photo(params.λ[1+floor(Int, mod(t, 364days)/day)], params)*min(μ_C, max(μ_NO₃, μ_NH₄))

        j_NH₄ = min(j̃_NH₄, μ*params.K_A*(N+params.N_struct))

        r = _r(T, μ, j_NO₃ + j_NH₄, params)

        dA = (μ - ν) * A / (60*60*24)
        dN = ((j_NO₃ + j_NH₄) / params.K_A - μ * (N + params.N_struct)) / (60*60*24)
        dC = ((p* (1 - e) - r) / params.K_A - μ * (C + params.C_struct)) / (60*60*24)

        A_new = A+dA*Δt 
        N_new = N+dN*Δt 
        C_new = C+dC*Δt

        if C_new < params.C_min
            A_new *= (1-(params.C_min - C)/params.C_struct)
            C_new = params.C_min
        end

        pp = (p-r)*A_new / (60*60*24*12*0.001) #gC/dm^2/hr to mmol C/s
        e *= p*A_new / (60*60*24*12*0.001)#mmol C/s
        ν *= params.K_A*A_new*(N_new + params.N_struct) / (60*60*24*14*0.001)#1/hr to mmol N/s
        j_NO₃ *= A_new / (60*60*24*14*0.001)#gN/dm^2/hr to mmol N/s
        j_NH₄ *= A_new / (60*60*24*14*0.001)#gN/dm^2/hr to mmol N/s

        du, dv, dw = 0.0, 0.0, 0.0
    else
        A_new, N_new, C_new, j_NO₃, j_NH₄, pp, e, ν, du, dv, dw = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end
    return (A = A_new, N = N_new, C = C_new, j_NO₃ = -j_NO₃, j_NH₄ = -j_NH₄, j_DIC = -pp, j_OXY=pp, e = e/6.56, ν = ν, u = du, v = dv, w = dw)
end

#fixed urel, T and S functions
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, irr::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, params.urel, params, Δt)
#fixed urel, T and S fields
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, params.urel, params, Δt)
#tracked u, T and S functions can not be done (like this) bcs same number of variables as main function
 #equations(x::AbstractFloat, y::AbstractFloat, z::Abstractparams.T(x, y, z, t)Float, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, v::AbstractFloat, w::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, sqrt(u^2+v^2+w^2), params, Δt)
#tracked u, T and S fields
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, v::AbstractFloat, w::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, T, S, irr, sqrt(u^2+v^2+w^2), params, Δt)

#parameters (some derived so need to be pre specified)
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

const defaults = (
    A_0 = 4.5,#6,# Growth rate adjustment parameter
    α = 4.15e-5 * 24 * 10^6 / (24 * 60 * 60),#3.75e-5 * 24 * 10^6 / (24 * 60 * 60),# photosynthetic efficiency 
    C_min = 0.01,# Minimal carbon reserve
    C_struct = 0.2,# Amount of carbon per unit dry weight of structural mass
    γ = 0.5,# Exudation parameter
    ϵ = 0.22,# Frond erosion parameter
    I_sat = 90 * 24*60*60/(10^6),#200 * 24 * 60 * 60 / (10^6),# Irradiance for maximal photosynthesis
    J_max = 1.4e-4 * 24,# Maximal nitrate uptake (gN/dm^2/h converted to gN/dm^2/day)
    K_A = .5,#0.6,# Structural dry weight per unit area
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

    j_NO₃_max = 10*0.5*24*14/(10^6), #Ahn 1998
    j_NH₄_max = 12*0.5*24*14/(10^6), #Ahn 1998

    uₐ = 0.72, #broch 2019
    uᵦ = 0.28,
    u₀ = 0.045,

    R_A=1.11e-4*24,
    R_B=5.57e-5*24
)

function gen_λ(lat)
    θ = 0.2163108 .+ 2 .* atan.(0.9671396 .* tan.(0.00860 .* ([1:365;] .- 186)))
    δ = asin.(0.39795 .* cos.(θ))
    p = 0.8333
    d = 24 .- (24 / pi) .* acos.((sin(p * pi / 180) .+ sin(lat * pi / 180) .* sin.(δ)) ./ (cos(lat * pi / 180) .* cos.(δ)))
    λ = diff(d)
    push!(λ, λ[end])
    return λ ./ findmax(λ)[1]
end

struct Particle
    #position
    x :: AbstractFloat
    y :: AbstractFloat
    z :: AbstractFloat
    #velocity
    u :: AbstractFloat
    v :: AbstractFloat
    w :: AbstractFloat
    #properties
    A :: AbstractFloat
    N :: AbstractFloat
    C :: AbstractFloat
    #feedback
    j_NO₃ :: AbstractFloat
    j_NH₄ :: AbstractFloat
    j_DIC :: AbstractFloat
    j_OXY :: AbstractFloat
    e :: AbstractFloat
    ν :: AbstractFloat
    #tracked fields
    NO₃  :: AbstractFloat
    NH₄  :: AbstractFloat
    PAR :: AbstractFloat
end

source_fields = (NO₃ = :NO₃, PAR=:PAR)

optional_source_fields = (NH₄ = :NH₄, )

#When Oceananigans PR in place going to simplify this specifcation
sink_fields = (NO₃ = (property=:j_NO₃, fallback=:N, fallback_scalefactor=(property=:A, constant=14*0.001/defaults.K_A)), )    
optional_sink_fields = (NH₄ = (property=:j_NH₄, fallback=:N, fallback_scalefactor=(property=:A, constant=14*0.001/defaults.K_A)), 
                                    DIC = (property=:j_DIC, fallback=:C, fallback_scalefactor=(property=:A, constant=14*0.001)),
                                    DOM = (property=:e, fallback=:A, fallback_scalefactor=0),#Rd_dom from LOBSTER, placeholder fallbacks because if these are taking away something has gone very wrong
                                    DD = (property=:ν, fallback=:A, fallback_scalefactor=0),
                                    OXY = (property=:j_OXY, fallback=:C, fallback_scalefactor=(property=:A, constant=14*0.001)))

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
    return StructArray{Particle}((x̄₀[1], x̄₀[2], x̄₀[3], zeros(n), zeros(n), zeros(n), x̄₀[4], x̄₀[5], x̄₀[6], [zeros(n) for i in 1:9]...))
end

@inline no_dynamics(args...) = nothing

function setup(n, x₀, y₀, z₀, A₀, N₀, C₀, latitude, density; T=nothing, S=nothing, urel=nothing, paramset=defaults, custom_dynamics=no_dynamics, optional_sources=(), optional_sinks=(), tracer_names=NamedTuple())
    if (!isnothing(T) && !isnothing(S) && isnothing(urel))
        throw(ArgumentError("T and S functions with tracked velocity fields not currently implimented"))
    elseif ((isnothing(T) && !isnothing(S)) | (!isnothing(T) && isnothing(S)))
        throw(ArgumentError("T and S must both be functions or both be tracked fields"))
    end

    particles = defineparticles((x₀=x₀, y₀=y₀, z₀=z₀, A₀=A₀, N₀=N₀, C₀=C₀), n)
    property_dependencies = (:A, :N, :C, :NO₃, :NH₄, :PAR)
    λ_arr=gen_λ(latitude)
    parameters = merge(paramset, (λ=λ_arr, ))
    tracked_properties = (:A, :N, :C, :j_NO₃, :j_NH₄, :j_DIC, :j_OXY, :e, :ν)

    if isnothing(T) property_dependencies = (property_dependencies..., :T, :S) else parameters = merge(parameters, (T=T, S=S)) end
    if isnothing(urel) property_dependencies = (property_dependencies..., :u, :v, :w) else parameters = merge(parameters, (urel = urel, )) end

    integral_properties = (:u, :v, :w) #for when I impliment dynamics

    sources = source_fields
    for tracer in optional_sources
        if tracer in keys(optional_source_fields)
            sources = merge(sources, NamedTuple{(tracer, )}((getproperty(optional_source_fields, tracer), )))
        else
            @warn "$tracer isn't an optional source field for SLatissima model"
        end
    end

    sinks = sink_fields
    for tracer in optional_sinks
        if tracer in keys(optional_sink_fields)
            sinks = merge(sinks, NamedTuple{(tracer, )}((getproperty(optional_sink_fields, tracer), )))
        else
            @warn "$tracer isn't an optional sink field for SLatissima model"
        end
    end

    sources = merge(sources, tracer_names)

    return Particles.setup(particles, equations, 
                                        property_dependencies,
                                        parameters,
                                        integral_properties, #changed all kelp model variables to diagnostic since it is much easier to enforce the extreme carbon limit correctly
                                        tracked_properties, 
                                        sources,
                                        sinks,
                                        density,
                                        custom_dynamics)
end
end

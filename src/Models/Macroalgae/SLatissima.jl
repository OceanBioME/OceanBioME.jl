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
using StructArrays, SugarKelp
using OceanBioME: Particles
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

@inline f_curr(u, params) = params.uₐ*(1-exp(-u/params.u₀))+params.uᵦ

function equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, params, Δt::AbstractFloat)
    if !iszero(A)
        irr/=3.99e-10*545e12/(1day) #W/m²/s to einstein/m²/day
        p=SugarKelp.p(T, irr, params)
        e=SugarKelp.e(C, params)
        ν = SugarKelp.ν(A, params)

        j_NO₃ = params.j_NO₃_max*f_curr(u, params)*((params.N_max-N)/(params.N_max-params.N_min))*NO₃/(params.k_NO₃+NO₃)
        j̃_NH₄ = params.j_NH₄_max*f_curr(u, params)*NH₄/(params.k_NH₄+NH₄)
        μ_NH₄ = j̃_NH₄/(params.K_A*(N+params.N_struct))
        μ_NO₃ = 1 - params.N_min/N
        μ_C = 1 - params.C_min/C

        μ = SugarKelp.f_area(A, params)*SugarKelp.f_temp(T, params)*SugarKelp.f_photo(params.λ[1+floor(Int, mod(t, 364days)/day)], params)*min(μ_C, max(μ_NO₃, μ_NH₄))

        j_NH₄ = min(j̃_NH₄, μ*params.K_A*(N+params.N_struct))

        r = SugarKelp.r(T, μ, j_NO₃ + j_NH₄, params.resp_model, params)

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
    return (A = A_new, N = N_new, C = C_new, j_NO₃ = j_NO₃, j_NH₄ = j_NH₄, pp = pp, e = e, ν = ν, u = du, v = dv, w = dw)
end

#fixed urel, T and S functions
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, irr::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, params.urel, params, Δt)
#fixed urel, T and S fields
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, params.urel, params, Δt)
#tracked u, T and S functions can not be done (like this) bcs same number of variables as main function
 #equations(x::AbstractFloat, y::AbstractFloat, z::Abstractparams.T(x, y, z, t)Float, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, v::AbstractFloat, w::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, params.T(x, y, z, t), params.S(x, y, z, t), irr, sqrt(u^2+v^2+w^2), params, Δt)
#tracked u, T and S fields
equations(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, A::AbstractFloat, N::AbstractFloat, C::AbstractFloat, NO₃::AbstractFloat, NH₄::AbstractFloat, T::AbstractFloat, S::AbstractFloat, irr::AbstractFloat, u::AbstractFloat, v::AbstractFloat, w::AbstractFloat, params, Δt::AbstractFloat) = equations(x, y, z, t, A, N, C, NO₃, NH₄, T, S, irr, sqrt(u^2+v^2+w^2), params, Δt)

const defaults = merge(SugarKelp.broch2013params, (
    j_NO₃_max = 10*0.5*24*14/(10^6), #Ahn 1998
    j_NH₄_max = 12*0.5*24*14/(10^6), #Ahn 1998
    k_NO₃ = SugarKelp.broch2013params.K_X,
    k_NH₄ = 1.3, #Fossberg 2018
    uₐ = 0.72, #broch 2019
    uᵦ = 0.28,
    u₀ = 0.045
))

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
    pp :: AbstractFloat
    e :: AbstractFloat
    ν :: AbstractFloat
    #tracked fields
    NO₃  :: AbstractFloat
    NH₄  :: AbstractFloat
    PAR :: AbstractFloat
end

source_fields = ((tracer=:NO₃, property=:NO₃, scalefactor=1.0), 
                        (tracer=:NH₄, property=:NH₄, scalefactor=1.0),
                        (tracer=:PAR, property=:PAR, scalefactor=1.0))
sink_fields = ((tracer=:NO₃, property=:j_NO₃, scalefactor=-1.0, fallback=:N, fallback_scalefactor=(property=:A, constant=14*0.001/defaults.K_A)), 
                    (tracer=:NH₄, property=:j_NH₄, scalefactor=-1.0, fallback=:N, fallback_scalefactor=(property=:A, constant=14*0.001/defaults.K_A)), 
                    (tracer=:DIC, property=:pp, scalefactor=-1.0, fallback=:C, fallback_scalefactor=(property=:A, constant=14*0.001)),
                    (tracer=:DOM, property=:e, scalefactor=1.0/6.56, fallback=:A, fallback_scalefactor=0),#Rd_dom from LOBSTER, placeholder fallbacks because if these are taking away something has gone very wrong
                    (tracer=:DD, property=:ν, scalefactor=1.0, fallback=:A, fallback_scalefactor=0))

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
    return StructArray{Particle}((x̄₀[1], x̄₀[2], x̄₀[3], zeros(n), zeros(n), zeros(n), x̄₀[4], x̄₀[5], x̄₀[6], [zeros(n) for i in 1:8]...))
end

@inline no_dynamics(args...) = nothing

function setup(n, x₀, y₀, z₀, A₀, N₀, C₀, latitude, density, T=nothing, S=nothing, urel=nothing, resp_model=2, paramset=defaults, custom_dynamics=no_dynamics, O₂=false)
    if (!isnothing(T) && !isnothing(S) && isnothing(urel))
        throw(ArgumentError("T and S functions with tracked velocity fields not currently implimented"))
    elseif ((isnothing(T) && !isnothing(S)) | (!isnothing(T) && isnothing(S)))
        throw(ArgumentError("T and S must both be functions or both be tracked fields"))
    end

    particles = defineparticles((x₀=x₀, y₀=y₀, z₀=z₀, A₀=A₀, N₀=N₀, C₀=C₀), n)
    property_dependencies = (:A, :N, :C, :NO₃, :NH₄, :PAR)
    λ_arr=SugarKelp.gen_λ(latitude)
    parameters = merge(paramset, (λ=λ_arr, resp_model=2))
    tracked_properties = (:A, :N, :C, :j_NO₃, :j_NH₄, :pp, :e, :ν)

    if isnothing(T) property_dependencies = (property_dependencies..., :T, :S) else parameters = merge(parameters, (T=T, S=S)) end
    if isnothing(urel) property_dependencies = (property_dependencies..., :u, :v, :w) else parameters = merge(parameters, (urel = urel, )) end

    integral_properties = (:u, :v, :w) #for when I impliment dynamics

    if O₂ merge(sink_fields, ((tracer=:OXY, property=:pp, scalefactor=1.0), )) end
    
    return Particles.setup(particles, equations, 
                                        property_dependencies,
                                        parameters,
                                        integral_properties, #changed all kelp model variables to diagnostic since it is much easier to enforce the extreme carbon limit correctly
                                        tracked_properties, 
                                        source_fields,
                                        sink_fields,
                                        density,
                                        custom_dynamics)
end
end

module Soetaert

using Roots, Oceananigans, KernelAbstractions
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Architectures: device
using Oceananigans: AdvectiveForcing
using Oceananigans.Forcings: maybe_constant_field
using Oceananigans.Advection: div_Uc

#These were determined by fitting data from a different model by the paper this is from
#Clearly they are not  closed (e.g. if anything is zero we get NaN or Inf)
#It is also very uncomfortable to log dimensional things + we get a nondimensional parameter p_whatever when we have "nondimensional"/dimensional value
#also slightly concerned that these are per day values but its unclear where I'd change the units

const defaults= (
    #https://aslopubs.onlinelibrary.wiley.com/doi/epdf/10.4319/lo.1996.41.8.1651
    λᵣᵣ = 2/year,# 1/year to 1/s
    λᵣ = 0.2/year,#s   
    Rdᵣᵣ = 0.1509,#mmol N/mmol C@inline Nᵣₑ_forcing(i, j, k, grid, clock, model_fields, params) = nothing
    Rdᵣ = 0.13,#mmol N/mmol C
    Rd_red = 6.56#106/16#mmol C/mmol N
)

@inline p_nit(Nₘ, Cₘ, O₂, NH₄, k) = min(1, exp(-1.9785+0.2261*log(Cₘ)*log(O₂)-0.0615*log(Cₘ)^2-0.0289*log(k)*log(NH₄)-0.36109*log(Cₘ)-0.0232*log(Cₘ)*log(NH₄))/Nₘ)
@inline p_denit(Cₘ, O₂, NO₃, k) = min(1, exp(-3.0790+1.7509*log(Cₘ)+0.0593*log(NO₃)^2-0.1923*log(Cₘ)^2+0.0604*log(k)^2+0.0662*log(O₂)*log(k))/Cₘ)
@inline p_anox(Cₘ, O₂, NO₃, k) = min(1, exp(-3.9476+2.6269*log(Cₘ)-0.2426*log(Cₘ)^2-1.3349*log(k)+0.1826*log(O₂)*log(k)-0.0143*log(NO₃)^2)/Cₘ)
#values too hight for d<100 so (for now) maxing at d=100
@inline p_soliddep(d) = 0.233*(982*max(d, 100)^(-1.548))^0.336

@inline _Nₘ(Nᵣᵣ, Nᵣ, params) = params.λᵣᵣ*Nᵣᵣ+params.λᵣ*Nᵣ
@inline _Cₘ(Nᵣᵣ, Nᵣ, params) = params.λᵣᵣ*params.Rdᵣᵣ*Nᵣᵣ+params.λᵣ*params.Rdᵣ*Nᵣ

@inline function mineralisation(Nᵣᵣ, Nᵣ, params)
    Nₘ=_Nₘ(Nᵣᵣ, Nᵣ, params)
    Cₘ=_Cₘ(Nᵣᵣ, Nᵣ, params)
    return Nₘ, Cₘ, Nₘ/(Nᵣᵣ+Nᵣ)
end

#params needs 2 fields: Nᵣᵣ and Nᵣ (like PAR in the lobster model), implimenting as discrete forcing so don't have to interpolate indices
@inline function sedimentNH₄(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(model_fields.Nᵣᵣ[i, j, 1],  model_fields.Nᵣ[i, j, 1], params)
    if (Nₘ>0 && Cₘ>0 && k > 0)
        return @inbounds Nₘ*(1-p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], k))
    else
        return 0
    end
end

@inline sedimentDIC(i, j, grid, clock, model_fields, params) = params.Rd_red*sedimentNH₄(i, j, grid, clock, model_fields, params)

@inline function sedimentNO₃(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(model_fields.Nᵣᵣ[i, j, 1],  model_fields.Nᵣ[i, j, 1], params)
    return @inbounds Nₘ*p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], k) - Cₘ*p_denit(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], k)*0.8
end

@inline function sedimentO₂(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(model_fields.Nᵣᵣ[i, j, 1],  model_fields.Nᵣ[i, j, 1], params)
    return @inbounds -(Cₘ*(1-p_denit(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], 1)-p_anox(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], 1)*p_soliddep(isa(params.d, Number) ? params.d : params.d[i, j])) + Nₘ*p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], 1)*2)
end

@inline POM_deposition(i, j, grid, clock, model_fields, params) = @inbounds - model_fields[params.POM][i, j, 1]*params.w*grid.Δzᵃᵃᶜ[1]

@inline eval_advection(i, j, k, grid, func) = func.advection_kernel_function(i, j, k, grid, func.advection_scheme, func.velocities, func.advected_field)
@inline getflux(i, j, k, grid, funcs) = length(funcs) < 11 ? ntuple(n->eval_advection(i, j, k, grid, funcs[n]), length(funcs)) : map(n -> eval_advection(i, j, k, grid, funcs[n]), 1:length(args))
@inline fPOM(i, j, k, grid, clock, model_fields, params) = - sum(getflux(i, j, k, grid, params.adv_funcs)).*params.f.*params.Δz

@inline Nᵣ_forcing(i, j, k, grid, clock, model_fields, params) = @inbounds - params.λᵣ*model_fields.Nᵣ[i, j, 1]
@inline Nᵣᵣ_forcing(i, j, k, grid, clock, model_fields, params) = @inbounds - params.λᵣᵣ*model_fields.Nᵣ[i, j, 1]

function setupsediment(grid, POM_fields::NamedTuple; POM_w = (DD = 200/day, D = 3.47e-5), parameters=defaults, Nᵣᵣᵢ=24.0, Nᵣᵢ=85.0, f_ref=0.1, f_fast=0.74, f_slow=0.26, carbonates=true)
    #fractions from Boudreau B. P. and Ruddick B. R. ( 1991) On a reactive continuum representation of organic matter diagenesis. Amer. J. Sci. 291, 507-538.
    #as used by https://reader.elsevier.com/reader/sd/pii/0016703796000130
    #additionally refractory fraction from https://reader.elsevier.com/reader/sd/pii/001670379500042X\
    #initial values from steady state from average deposition rate after 1 year in model without sediment (about 0.01mmolN/m³ for DD and 0 for D)
    if !(f_fast+f_slow==1)
        throw(ArgumentError("Reactivity fractions do not sum to 1"))
    end
    Nᵣᵣ=Oceananigans.Fields.Field{Center, Center, Center}(grid; indices=(:, :, 1:1))
    Nᵣ=Oceananigans.Fields.Field{Center, Center, Center}(grid; indices=(:, :, 1:1))
    Nᵣₑ=Oceananigans.Fields.Field{Center, Center, Center}(grid; indices=(:, :, 1:1))

    #would like to user smaller grid to concerve memeory but breaks output writer
    Nᵣᵣ .= Nᵣᵣᵢ; Nᵣ .= Nᵣᵢ; Nᵣₑ .= 0.0
    d = grid.zᵃᵃᶠ[1]
    parameters = merge(parameters, (; d, Nᵣᵣ, Nᵣ, Nᵣₑ, f_ref, f_fast, f_slow, POM_fields, POM_w))

    if carbonates DIC_bc = (DIC = FluxBoundaryCondition(sedimentDIC, discrete_form=true, parameters=parameters), ) else DIC_bc = () end
    
    fPOM_forcings=[]
    for (POM_name, POM) in pairs(POM_fields)
        u, v, w = maybe_constant_field.((0.0, 0.0, POM_w[POM_name]))
        velocities = (; u, v, w)
        push!(fPOM_forcings, AdvectiveForcing(velocities, WENO(;grid), div_Uc, POM))
    end

    return (
        boundary_conditions = merge((
                NO₃ = FluxBoundaryCondition(sedimentNO₃, discrete_form=true, parameters=parameters),
                NH₄ = FluxBoundaryCondition(sedimentNH₄, discrete_form=true, parameters=parameters),
                OXY = FluxBoundaryCondition(sedimentO₂, discrete_form=true, parameters=parameters)
            ), 
            DIC_bc
        ),
        forcing = (
            Nᵣ = (Forcing(Nᵣ_forcing, discrete_form=true, parameters=parameters), Forcing(fPOM, discrete_form=true, parameters=(adv_funcs = fPOM_forcings, f = (1 - f_ref)*f_slow, Δz = grid.Δzᵃᵃᶜ[1]))), 
            Nᵣᵣ = (Forcing(Nᵣᵣ_forcing, discrete_form=true, parameters=parameters), Forcing(fPOM, discrete_form=true, parameters=(adv_funcs = fPOM_forcings, f = (1 - f_ref)*f_fast, Δz = grid.Δzᵃᵃᶜ[1]))), 
            Nᵣₑ = Forcing(fPOM, discrete_form=true, parameters=(adv_funcs = fPOM_forcings, f = f_ref, Δz = grid.Δzᵃᵃᶜ[1]))
        ),
        auxiliary_fields = (Nᵣᵣ=Nᵣᵣ, Nᵣ=Nᵣ, Nᵣₑ=Nᵣₑ),
        parameters = parameters
    )
end
end
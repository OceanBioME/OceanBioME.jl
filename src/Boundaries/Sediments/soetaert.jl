"""
Soetaert

Soetaert et al. 2000 vertically integrated sediment diagenesis model.

XXX Currently does not work properly and will result in nitrogen leak!

Notes
=======
- In the paper `p` functions are determined from fitting data and are not closed (e.g. if anything is zero they go to NaN or Inf)
- Additionally theres a lot of `log` and `exp` of dimensional quantities 
- Units may also be wrong
"""
module Soetaert

using Roots, Oceananigans, KernelAbstractions
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Architectures: device
using Oceananigans: AdvectiveForcing
using Oceananigans.Forcings: maybe_constant_field
using Oceananigans.Advection: div_Uc

#####
##### Default parameter definition
#####

const defaults= (
    #https://aslopubs.onlinelibrary.wiley.com/doi/epdf/10.4319/lo.1996.41.8.1651
    λᵣᵣ = 2*day/year,# 1/year to 1/s
    λᵣ = 0.2*day/year,#s   
    Rdᵣᵣ = 0.1509,#mmol N/mmol C@inline Nᵣₑ_forcing(i, j, k, grid, clock, model_fields, params) = nothing
    Rdᵣ = 0.13,#mmol N/mmol C
    Rd_red = 6.56#106/16#mmol C/mmol N
)

#####
##### Fractional nitrifcation/denitrifcation/anoxic mineralisation and mineralisation rate calculations
#####
@inline p_nit(Nₘ, Cₘ, O₂, NH₄, k) = exp(-1.9785+0.2261*log(Cₘ)*log(O₂)-0.0615*log(Cₘ)^2-0.0289*log(k)*log(NH₄)-0.36109*log(Cₘ)-0.0232*log(Cₘ)*log(NH₄))/Nₘ

@inline p_denit(Cₘ, O₂, NO₃, k) = exp(-3.0790+1.7509*log(Cₘ)+0.0593*log(NO₃)^2-0.1923*log(Cₘ)^2+0.0604*log(k)^2+0.0662*log(O₂)*log(k))/Cₘ

@inline p_anox(Cₘ, O₂, NO₃, k) = exp(-3.9476+2.6269*log(Cₘ)-0.2426*log(Cₘ)^2-1.3349*log(k)+0.1826*log(O₂)*log(k)-0.0143*log(NO₃)^2)/Cₘ

# XXX values too hight for d<100 so (for now) maxing at d=100
@inline p_soliddep(d) = 0.233*(982*max(d, 100)^(-1.548))^0.336

@inline _Nₘ(Nᵣᵣ, Nᵣ, params) = params.λᵣᵣ*Nᵣᵣ+params.λᵣ*Nᵣ #mmol N/day
@inline _Cₘ(Nᵣᵣ, Nᵣ, params) = params.λᵣᵣ*params.Rdᵣᵣ*Nᵣᵣ+params.λᵣ*params.Rdᵣ*Nᵣ #mmol C/day

@inline function mineralisation(Nᵣᵣ, Nᵣ, params)
    Nₘ=_Nₘ(Nᵣᵣ, Nᵣ, params)
    Cₘ=_Cₘ(Nᵣᵣ, Nᵣ, params)
    return Nₘ, Cₘ, Nₘ/(Nᵣᵣ+Nᵣ) #mmol N/day, mmol C/day, 1/day
end

#####
##### Sediment boundary condition functions
#####

@inline function sedimentNH₄(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(model_fields.Nᵣᵣ[i, j, 1],  model_fields.Nᵣ[i, j, 1], params)
    return @inbounds (Nₘ*(1-p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], k)) + Cₘ*p_denit(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], k)*0.8)/day
end

@inline sedimentDIC(i, j, grid, clock, model_fields, params) = params.Rd_red*sedimentNH₄(i, j, grid, clock, model_fields, params)

@inline function sedimentNO₃(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(model_fields.Nᵣᵣ[i, j, 1],  model_fields.Nᵣ[i, j, 1], params)
    return @inbounds (Nₘ*p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], k) - Cₘ*p_denit(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], k)*0.8)/day
end

@inline function sedimentO₂(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(model_fields.Nᵣᵣ[i, j, 1],  model_fields.Nᵣ[i, j, 1], params)
    return @inbounds -(Cₘ*(1-p_denit(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], 1)-p_anox(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], 1)*p_soliddep(isa(params.d, Number) ? params.d : params.d[i, j])) + Nₘ*p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], 1)*2)/day
end

#####
##### Sediment prognostics forcing 
#####
@inline eval_advection(i, j, k, grid, func) = func.advection_kernel_function(i, j, k, grid, func.advection_scheme, func.velocities, func.advected_field)
@inline getflux(i, j, k, grid, funcs) = length(funcs) < 11 ? ntuple(n->eval_advection(i, j, k, grid, funcs[n]), length(funcs)) : map(n -> eval_advection(i, j, k, grid, funcs[n]), 1:length(args))
@inline fPOM(i, j, k, grid, clock, model_fields, params) = - sum(getflux(i, j, k, grid, params.adv_funcs)).*params.f.*params.Δz

@inline Nᵣ_forcing(i, j, k, grid, clock, model_fields, params) = @inbounds getflux(i, j, grid, params.adv_funcs).*params.f.*params.Δz - params.λᵣ*model_fields.Nᵣ[i, j, 1]/day
@inline Nᵣᵣ_forcing(i, j, k, grid, clock, model_fields, params) = @inbounds getflux(i, j, grid, params.adv_funcs).*params.f.*params.Δz - params.λᵣᵣ*model_fields.Nᵣ[i, j, 1]/day

"""
    Boundaries.Sediments.Soetaert.setup(grid, POM_fields::NamedTuple; 
                                                                POM_w = (DD = 200/day, D = 3.47e-5), 
                                                                parameters=defaults, 
                                                                Nᵣᵣᵢ=24.0, 
                                                                Nᵣᵢ=85.0, 
                                                                f_ref=0.1, 
                                                                f_fast=0.74,
                                                                f_slow=0.26, 
                                                                carbonates=true)

Returns the boundary conditions, forcing, auxiliary fields, and parameters for the [Soetaert2000](@cite)
vertically integrated sediment diagenesis model. 

Arguments
=========
* `grid`: Oceananigans grid.
* `POM_fields`: NamedTuple of Oceananigans fields which contribute organic matter to the sediment
Keyword arguments
=================
* `POM_w`: NamedTuple of the sinking speed (at the bottom of the domain) for each sinking tracer
* `Nᵣᵣᵢ` and `Nᵣᵢ`: initial values for the fast and slow degrading sediment components
* `f_ref`, `f_fast`, and `f_slow`: fraction of sedimened organic matter which is: refractory, fast degrading and slow degrading. 
    `f_ref` is a fraction of the total depositional flux where as the others are a fraction of the remaining flux
* `carbonates`: as a temporary measure this is used to determine if a DIC boundary condition should be provided

Notes
=======
Fast/slow/refractory fractions from [Boudreau1991](@cite) and [Tromp1995](@cite) as used by [Soetaert2000](@cite) 
and initial values from steadyt state average deposition after 1 year without sediment 
in sub polar region (approx 0.01 mmolN/m³ for DD and 0 for D).
"""
function setup(grid, POM_fields::NamedTuple; 
                                      POM_w = (DD = 200/day, D = 3.47e-5), 
                                      parameters=defaults, 
                                      Nᵣᵣᵢ=24.0, 
                                      Nᵣᵢ=85.0, 
                                      f_ref=0.1, 
                                      f_fast=0.74,
                                      f_slow=0.26, 
                                      carbonates=true)

    if !(f_fast+f_slow==1)
        throw(ArgumentError("Reactivity fractions do not sum to 1"))
    end
    Nᵣᵣ=Oceananigans.Fields.Field{Center, Center, Center}(grid; indices=(:, :, 1:1))
    Nᵣ=Oceananigans.Fields.Field{Center, Center, Center}(grid; indices=(:, :, 1:1))
    Nᵣₑ=Oceananigans.Fields.Field{Center, Center, Center}(grid; indices=(:, :, 1:1))

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
            Nᵣ = Forcing(Nᵣ_forcing, discrete_form=true, parameters=merge(parameters,(adv_funcs = fPOM_forcings, f = (1 - f_ref)*f_slow, Δz = grid.Δzᵃᵃᶜ[1]))), 
            Nᵣᵣ = Forcing(Nᵣᵣ_forcing, discrete_form=true, parameters=merge(parameters, (adv_funcs = fPOM_forcings, f = (1 - f_ref)*f_fast, Δz = grid.Δzᵃᵃᶜ[1]))), 
            Nᵣₑ = Forcing(fPOM, discrete_form=true, parameters=(adv_funcs = fPOM_forcings, f = f_ref, Δz = grid.Δzᵃᵃᶜ[1]))
        ),
        auxiliary_fields = (Nᵣᵣ=Nᵣᵣ, Nᵣ=Nᵣ, Nᵣₑ=Nᵣₑ),
        parameters = parameters
    )
end
end
module Soetaert

using Roots, Oceananigans, KernelAbstractions
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years
using Oceananigans.Architectures: device

#These were determined by fitting data from a different model by the paper this is from
#Clearly they are not  closed (e.g. if anything is zero we get NaN or Inf)
#It is also very uncomfortable to log dimensional things + we get a nondimensional parameter p_whatever when we have "nondimensional"/dimensional value
#also slightly concerned that these are per day values but its unclear where I'd change the units

const defaults= (
    #https://aslopubs.onlinelibrary.wiley.com/doi/epdf/10.4319/lo.1996.41.8.1651
    λᵣᵣ = 2/year,# 1/year to 1/s
    λᵣ = 0.2/year,#s   
    Rdᵣᵣ = 0.1509,#mmol N/mmol C
    Rdᵣ = 0.13,#mmol N/mmol C
    Rd_red = 106/16#mmol C/mmol N
)

@inline p_nit(Nₘ, Cₘ, O₂, NH₄, k) = min(1, exp(-1.9785+0.2261*log(Cₘ)*log(O₂)-0.0615*log(Cₘ)^2-0.0289*log(k)*log(NH₄)-0.36109*log(Cₘ)-0.0232*log(Cₘ)*log(NH₄))/Nₘ)
@inline p_denit(Cₘ, O₂, NO₃, k) = min(1, exp(-3.0790+1.7509*log(Cₘ)+0.0593*log(NO₃)^2-0.1923*log(Cₘ)^2+0.0604*log(k)^2+0.0662*log(O₂)*log(k))/Cₘ)
@inline p_anox(Cₘ, O₂, NO₃, k) = min(1, exp(-3.9476+2.6269*log(Cₘ)-0.2426*log(Cₘ)^2-1.3349*log(k)+0.1826*log(O₂)*log(k)-0.0143*log(NO₃)^2)/Cₘ)
#values too hight for d<100 so (for now) maxing at d=100
@inline p_soliddep(d) = 0.233*(982*max(d, 100)^(-1.548))^0.336

@inline _Nₘ(Nᵣᵣ, Nᵣ, params) = params.λᵣᵣ*Nᵣᵣ+params.λᵣ*Nᵣ
@inline _Cₘ(Nᵣᵣ, Nᵣ, params) = params.λᵣᵣ*params.Rdᵣᵣ*Nᵣᵣ+params.λᵣ*params.Rdᵣ*Nᵣ

@inline function mineralisation(i, j, params)
    Nᵣᵣ=params.Nᵣᵣ[i, j, 1]
    Nᵣ=params.Nᵣ[i, j, 1]
    Nₘ=_Nₘ(Nᵣᵣ, Nᵣ, params)
    Cₘ=_Cₘ(Nᵣᵣ, Nᵣ, params)
    return Nₘ, Cₘ, Nₘ/(Nᵣᵣ+Nᵣ)
end

#params needs 2 fields: Nᵣᵣ and Nᵣ (like PAR in the lobster model), implimenting as discrete forcing so don't have to interpolate indices
@inline function sedimentNH₄(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(i, j, params)
    if (Nₘ>0 && Cₘ>0 && k > 0)
        return @inbounds Nₘ*(1-p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], k))#max just incase the odd p_nit etc give unphysical values (this can not be a flux out of the system as it is irreversible chemistry)
    else
        return 0
    end
end

@inline sedimentDIC(i, j, grid, clock, model_fields, params) = params.Rd_red*sedimentNH₄(i, j, grid, clock, model_fields, params)

@inline function sedimentNO₃(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(i, j, params)
    if (Nₘ>0 && Cₘ>0 && k > 0)
        return @inbounds Nₘ*p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], k) - Cₘ*p_denit(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], k)*0.8
    else
        return 0
    end
end

@inline function sedimentO₂(i, j, grid, clock, model_fields, params)
    Nₘ, Cₘ, k = mineralisation(i, j, params)
    if (Nₘ>0 && Cₘ>0 && k > 0)
        return @inbounds -(Cₘ*(1-p_denit(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], 1)-p_anox(Cₘ, model_fields.OXY[i, j, 1], model_fields.NO₃[i, j, 1], 1)*p_soliddep(isa(params.d, Number) ? params.d : params.d[i,j])) + Nₘ*p_nit(Nₘ, Cₘ, model_fields.OXY[i, j, 1], model_fields.NH₄[i, j, 1], 1)*2)
    else
        return 0
    end
end

@inline sedimentD(i, j, grid, clock, model_fields, params) = @inbounds -model_fields.D[i, j, 1]*params.w_slow*tanh(max(-params.d/params.λ, 0))
@inline sedimentDD(i, j, grid, clock, model_fields, params) = @inbounds -model_fields.DD[i, j, 1]*params.w_fast*tanh(max(-params.d/params.λ, 0))

@kernel function _integrate_sediment!(D, DD, Nᵣᵣ, Nᵣ, Nᵣₑ, Δt, params)
    i, j = @index(Global, NTuple) 
    N_dep = (params.w_slow*tanh(max(-params.d/params.λ, 0)) * D[i, j, 1] + params.w_fast*tanh(max(-params.d/params.λ, 0)) * DD[i, j, 1])*2#no idea why we need this factor of 2 but nitrogen budget doesnt add up otherwise
    Nᵣᵣ[i, j, 1] += (N_dep*params.f_fast*(1-params.f_ref) - params.λᵣᵣ*Nᵣᵣ[i, j, 1]) *Δt
    Nᵣ[i, j, 1] += (N_dep*params.f_slow*(1-params.f_ref) - params.λᵣ*Nᵣ[i, j, 1]) *Δt
    Nᵣₑ[i, j, 1] += N_dep*params.f_ref*Δt
end

function integrate_sediment!(sim, params)
    calculation = Oceananigans.Utils.launch!(sim.model.architecture, sim.model.grid, :xy, _integrate_sediment!,
                                   sim.model.tracers.D, sim.model.tracers.DD, sim.model.auxiliary_fields.Nᵣᵣ, sim.model.auxiliary_fields.Nᵣ, sim.model.auxiliary_fields.Nᵣₑ, sim.Δt, params,
                                   dependencies = Event(device(sim.model.architecture)))
    wait(device(sim.model.architecture), calculation)
end

function setupsediment(grid; w_fast=200/day, w_slow=3.47e-5, λ=1.0, parameters=defaults, Nᵣᵣᵢ=24.0, Nᵣᵢ=85.0, f_ref=0.1, f_fast=0.74, f_slow=0.26, carbonates=true)
    @warn "Sediment model will corrently produce incorrect results with non-zero water velocity (currently the deposition simply the concentration*sinking speed which is incorrect if detritus has other w components"
    #fractions from Boudreau B. P. and Ruddick B. R. ( 1991) On a reactive continuum representation of organic matter diagenesis. Amer. J. Sci. 291, 507-538.
    #as used by https://reader.elsevier.com/reader/sd/pii/0016703796000130
    #additionally refractory fraction from https://reader.elsevier.com/reader/sd/pii/001670379500042X
    #don't actually need to replicate the size just nuimber of nodes but doing incase we interpolate in the future
    #bottomgrid=RectilinearGrid(size=(grid.Nx, grid.Ny, 1), extent=(grid.Lx, grid.Ly, 1.0))
    #initial values from steady state from average deposition rate after 1 year in model without sediment (about 0.01mmolN/m³ for DD and 0 for D)
    if !(f_fast+f_slow==1)
        throw(ArgumentError("Reactivity fractions do not sum to 1"))
    end
    Nᵣᵣ=Oceananigans.Fields.Field{Center, Center, Center}(grid)#bottomgrid)
    Nᵣ=Oceananigans.Fields.Field{Center, Center, Center}(grid)#bottomgrid)
    Nᵣₑ=Oceananigans.Fields.Field{Center, Center, Center}(grid)#bottomgrid)

    #would like to user smaller grid to concerve memeory but breaks output writer
    Nᵣᵣ .= Nᵣᵣᵢ; Nᵣ .= Nᵣᵢ; Nᵣₑ .= 0.0
    parameters = merge(parameters, (Nᵣᵣ = Nᵣᵣ, Nᵣ = Nᵣ, Nᵣₑ=Nᵣₑ, d = grid.zᵃᵃᶠ[1], λ=λ, f_ref=f_ref, f_fast=f_fast, f_slow=f_slow, w_fast = w_fast, w_slow = w_slow))

    if carbonates DIC_bc = (DIC = FluxBoundaryCondition(sedimentDIC, discrete_form=true, parameters=parameters), ) else DIC_bc = () end

    return (
        boundary_conditions = merge((
            NO₃ = FluxBoundaryCondition(sedimentNO₃, discrete_form=true, parameters=parameters),
            NH₄ = FluxBoundaryCondition(sedimentNH₄, discrete_form=true, parameters=parameters),
            OXY = FluxBoundaryCondition(sedimentO₂, discrete_form=true, parameters=parameters),
            D = FluxBoundaryCondition(sedimentD, discrete_form=true, parameters=parameters),
            DD = FluxBoundaryCondition(sedimentDD, discrete_form=true, parameters=parameters)
        ), DIC_bc),
        auxiliary_fields = (Nᵣᵣ=Nᵣᵣ, Nᵣ=Nᵣ, Nᵣₑ=Nᵣₑ),
        callback =  Callback(integrate_sediment!, IterationInterval(1), parameters),
        parameters = parameters
    )
end
end
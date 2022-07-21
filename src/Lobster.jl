module Lobster
#https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2004JC002588
using Oceananigans

include("parameters/lobster.jl")
import BGC: BGCModel, AirSeaFlux

@inline p(P, D, params) = params.p̃*P/(params.p̃*P+(1-params.p̃)*D)
@inline G_d(P, Z, D, params) = params.g_z*(1-p(P, D, params))*D*Z/(params.K_z+P*p(P, D, params)+(1-p(P, D, params))*D)
@inline G_p(P, Z, D, params) = params.g_z*p(P, D, params)*P*Z/(params.K_z+P*p(P, D, params)+(1-p(P, D, params))*D)

#Limiting equations
@inline L_I(PAR, params) = 1-exp(-PAR/params.K_par) # reference Levy 2001
@inline L_NO₃(NO₃, NH₄, params) = NO₃*exp(-params.ψ*NH₄)/(NO₃+params.K_no₃)
@inline L_NH₄(NH₄, params) = max(NH₄/(NH₄+params.K_nh₄),0)    # NH₄/(NH₄+K_nh₄)

#source functions
Z_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = params.a_z*(G_d(P, Z, D, params) + G_p(P, Z, D, params)) - params.m_z*Z^2 - params.μ_z*Z
D_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = (1-params.f_d)*(1-params.a_z)*(G_d(P, Z, D, params)+G_p(P, Z, D, params)) + (1-params.f_d)*params.m_p*P^2 - G_d(P, Z, D, params) + params.f_z*params.m_z*Z^2 - params.μ_d*D #- aggreg_D2DD(z, D, DD) + aggreg_DOM2D(z, D, DOM) 
DD_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = params.f_d*(1-params.a_z)*(G_d(P, Z, D, params)+G_p(P, Z, D, params)) + params.f_d*params.m_p*P^2 + (1-params.f_z)*params.m_z*Z^2 - params.μ_dd*DD #+ aggreg_D2DD(z, D, DD) + aggreg_DOM2DD(z, DD, DOM) 
P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = (1-params.γ)*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P - G_p(P, Z, D, params) - params.m_p*P^2
NO₃_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = -params.μ_p*L_I(PAR, params)*L_NO₃(NO₃, NH₄, params)*P + params.μ_n*NH₄ #+ delta(t-t_inject)*exp(-(z+50)^2/10^2)*NO₃_flux
NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = params.α_p*params.γ*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P - params.μ_p*L_I(PAR, params)*L_NH₄(NH₄, params)*P - params.μ_n*NH₄ + params.α_z*params.μ_z*Z + params.α_d*params.μ_d*D + params.α_dd*params.μ_dd*DD + params.μ_dom*DOM + (1-params.Rd_phy/params.Rd_dom)*((1-params.α_p)*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P+(1-params.α_z)*params.μ_z*Z+(1-params.α_d)*params.μ_d*D+(1-params.α_dd)*params.μ_dd*DD) #+ exp(-(z+50)^2/10^2)*delta(t-t_inject)*NH₄_flux
DOM_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = (1-params.α_p)*params.γ*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P - params.μ_dom*DOM +(1-params.α_z)*params.μ_z*Z+(1-params.α_d)*params.μ_d*D +(1-params.α_dd)*params.μ_dd*DD - (1-params.Rd_phy/params.Rd_dom)*((1-params.α_p)*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P+(1-params.α_z)*params.μ_z*Z+(1-params.α_d)*params.μ_d*D+(1-params.α_dd)*params.μ_dd*DD) #- aggreg_DOM2D(z, D, DOM)- aggreg_DOM2DD(z, DD, DOM)
DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = -params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*params.Rd_phy*(1+params.ρ_caco3)*P+params.α_p*params.γ*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*params.Rd_phy*P+params.α_z*params.μ_z*params.Rd_phy*Z+params.α_d*params.μ_d*params.Rd_phy*D+params.α_dd*params.μ_dd*params.Rd_phy*DD+params.μ_dom*DOM*params.Rd_dom #+ delta(z-grid.zᵃᵃᶜ[Nz])*air_sea_flux(16, 2.002e-3, 2.311e-3, 8.0)/grid.Δzᵃᵃᶜ
ALK_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = params.μ_p*L_I(PAR, params)*L_NO₃(NO₃, NH₄, params)*P-2*params.ρ_caco3*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*params.Rd_phy*P

getPAR(PAR::Field, x, y, z, t) = Oceananigans.Fields.interpolate(PAR, Center(), Center(), Center(), PAR.grid, x, y, z)
getPAR(PAR, x, y, z, t) = PAR(x, y, z, t)

Z_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = Z_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
D_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = D_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
DD_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = DD_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) =P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
NO₃_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = NO₃_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
DOM_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = DOM_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
ALK_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = ALK_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)

function setup(grid, parameters, forcings=(T=nothing, S=nothing, PAR=nothing), fluxboundaries=(dic=AirSeaFlux.dic, ))
    if keys(forcings) != (:T, :S, :PAR)
        throw(ArgumentError("Need to provide all external forcings (T, S, PAR), either pass 'nothing' for a field to be used, a function of x, y, z, t, or an externally defined field which must also be an auxililiary field of the model. Currently T and S must be the same types."))
    end
    #setup forcings
    slip_vel_D = Oceananigans.Architectures.arch_array(grid.architecture, zeros(0:grid.Nx+2,0:grid.Ny+2,0:grid.Nz+2))
    slip_vel_DD = Oceananigans.Architectures.arch_array(grid.architecture, zeros(0:grid.Nx+2,0:grid.Ny+2,0:grid.Nz+2))
    @simd for k=0:grid.Nz-2           #0:Nz+2   3:Nz-2
        #slip_vel_D[:,:,k].+=V_d*(tanh(max(-grid.zᵃᵃᶠ[k]/λ,0)))#*tanh(max((grid.zᵃᵃᶠ[k]+Lz)/λ,0)));or pass 
        @inbounds slip_vel_D[:,:,k].+=parameters.V_d*tanh(max(-grid.zᵃᵃᶠ[k]/parameters.λ,0))#*tanh(max((grid.zᵃᵃᶠ[k]+Lz)/λ,0));
        @inbounds slip_vel_DD[:,:,k].+=parameters.V_dd*tanh(max(-grid.zᵃᵃᶠ[k]/parameters.λ,0))#*tanh(max((grid.zᵃᵃᶠ[k]+Lz)/λ,0));
    end
    D_RHS_slip = AdvectiveForcing(WENO5(; grid), w=slip_vel_D)
    DD_RHS_slip = AdvectiveForcing(WENO5(; grid), w=slip_vel_DD)

    if isnothing(forcings.PAR)
        throw(ArgumentError("PAR should not be a tracer field, pass fuinction of x, y, z, t, or an auxiliary field."))
    end

    function_field_dependencies = (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM, :DIC, :ALK)
    function_parameters = merge(parameters, (PAR = forcings.PAR, ))#this can be an auxiliary field or a function of x, y, z, t

    Z_RHS = Forcing(Z_forcing, field_dependencies = function_field_dependencies, parameters = function_parameters)
    D_RHS_nonslip = Forcing(D_forcing_nonslip, field_dependencies = function_field_dependencies, parameters = function_parameters)
    DD_RHS_nonslip = Forcing(DD_forcing_nonslip, field_dependencies = function_field_dependencies, parameters = function_parameters)
    P_RHS = Forcing(P_forcing, field_dependencies = function_field_dependencies, parameters = function_parameters)
    NO₃_RHS = Forcing(NO₃_forcing, field_dependencies = function_field_dependencies, parameters = function_parameters)
    NH₄_RHS = Forcing(NH₄_forcing, field_dependencies = function_field_dependencies, parameters = function_parameters)
    DOM_RHS = Forcing(DOM_forcing, field_dependencies = function_field_dependencies, parameters = function_parameters)
    DIC_RHS = Forcing(DIC_forcing, field_dependencies = function_field_dependencies, parameters = function_parameters)
    ALK_RHS = Forcing(ALK_forcing, field_dependencies = function_field_dependencies, parameters = function_parameters)

    #setup boundary conditions
    NO₃_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))
    NH₄_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))
    P_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))  
                                
    Z_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))
    D_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))
    DD_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))                                
    DOM_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))
    
    if (isnothing(forcings.T) & isnothing(forcings.S))
        bcs_field_dependencies = (:DIC, :ALK, :T, :S)
        bcs_parameters = parameters
    elseif isnothing(forcings.T)
        bcs_field_dependencies = (:DIC, :ALK, :T)
        bcs_parameters = merge(parameters, (S=forcings.S, ))
    elseif isnothing(forcings.S)
        bcs_field_dependencies = (:DIC, :ALK, :S)
        bcs_parameters = merge(parameters, (T=forcings.T, ))
    else
        bcs_field_dependencies = (:DIC, :ALK)
        bcs_parameters = merge(parameters, (T=forcings.T, S=forcings.S))
    end

    dic_bc = FluxBoundaryCondition(fluxboundaries.dic, field_dependencies = bcs_field_dependencies, parameters = bcs_parameters)
    DIC_bcs = FieldBoundaryConditions(top = dic_bc, bottom = FluxBoundaryCondition(0))             #airseaflux_bc, #FluxBoundaryCondition(-1.3e-4)                                
    ALK_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))

   return BGCModel((:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM, :DIC, :ALK), #tracers
            (Z=Z_RHS, D=(D_RHS_nonslip, D_RHS_slip), DD=(DD_RHS_nonslip, DD_RHS_slip), P=P_RHS, NO₃=NO₃_RHS, NH₄=NH₄_RHS, DOM=DOM_RHS, DIC=DIC_RHS, ALK=ALK_RHS), #forcing
            (P=P_bcs, Z=Z_bcs, D=D_bcs, DD=DD_bcs, NO₃=NO₃_bcs, NH₄=NH₄_bcs, DOM=DOM_bcs, DIC=DIC_bcs, ALK=ALK_bcs)) #boundaries
end
end # module

module LOBSTER
#https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2004JC002588
using Oceananigans
include("parameters/lobster.jl")

@inline p(P, D, params) = params.p̃*P/(params.p̃*P+(1-params.p̃)*D)
@inline G_d(P, Z, D, params) = params.g_z*(1-p(P, D, params))*D*Z/(params.K_z+P*p(P, D, params)+(1-p(P, D, params))*D)
@inline G_p(P, Z, D, params) = params.g_z*p(P, D, params)*P*Z/(params.K_z+P*p(P, D, params)+(1-p(P, D, params))*D)

#Limiting equations
@inline L_I(PAR, params) = 1-exp(-PAR/params.K_par) # reference Levy 2001
@inline L_NO₃(NO₃, NH₄, params) = NO₃*exp(-params.ψ*NH₄)/(NO₃+params.K_no₃)
@inline L_NH₄(NH₄, params) = max(NH₄/(NH₄+params.K_nh₄),0)    # NH₄/(NH₄+K_nh₄)

#source functions
Z_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR::AbstractFloat, params) = params.a_z*(G_d(P, Z, D, params) + G_p(P, Z, D, params)) - params.m_z*Z^2 - params.μ_z*Z
D_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR::AbstractFloat, params) = (1-params.f_d)*(1-params.a_z)*(G_d(P, Z, D, params)+G_p(P, Z, D, params)) + (1-params.f_d)*params.m_p*P^2 - G_d(P, Z, D, params) + params.f_z*params.m_z*Z^2 - params.μ_d*D #- aggreg_D2DD(z, D, DD) + aggreg_DOM2D(z, D, DOM) 
DD_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR::AbstractFloat, params) = params.f_d*(1-params.a_z)*(G_d(P, Z, D, params)+G_p(P, Z, D, params)) + params.f_d*params.m_p*P^2 + (1-params.f_z)*params.m_z*Z^2 - params.μ_dd*DD #+ aggreg_D2DD(z, D, DD) + aggreg_DOM2DD(z, DD, DOM) 
P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR::AbstractFloat, params) = (1-params.γ)*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P - G_p(P, Z, D, params) - params.m_p*P^2
NO₃_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR::AbstractFloat, params) = -params.μ_p*L_I(PAR, params)*L_NO₃(NO₃, NH₄, params)*P + params.μ_n*NH₄ #+ delta(t-t_inject)*exp(-(z+50)^2/10^2)*NO₃_flux
NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR::AbstractFloat, params) = params.α_p*params.γ*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P - params.μ_p*L_I(PAR, params)*L_NH₄(NH₄, params)*P - params.μ_n*NH₄ + params.α_z*params.μ_z*Z + params.α_d*params.μ_d*D + params.α_dd*params.μ_dd*DD + params.μ_dom*DOM + (1-params.Rd_phy/params.Rd_dom)*((1-params.α_p)*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P+(1-params.α_z)*params.μ_z*Z+(1-params.α_d)*params.μ_d*D+(1-params.α_dd)*params.μ_dd*DD) #+ exp(-(z+50)^2/10^2)*delta(t-t_inject)*NH₄_flux
DOM_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR::AbstractFloat, params) = (1-params.α_p)*params.γ*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P - params.μ_dom*DOM +(1-params.α_z)*params.μ_z*Z+(1-params.α_d)*params.μ_d*D +(1-params.α_dd)*params.μ_dd*DD - (1-params.Rd_phy/params.Rd_dom)*((1-params.α_p)*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*P+(1-params.α_z)*params.μ_z*Z+(1-params.α_d)*params.μ_d*D+(1-params.α_dd)*params.μ_dd*DD) #- aggreg_DOM2D(z, D, DOM)- aggreg_DOM2DD(z, DD, DOM)

DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = -params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*params.Rd_phy*(1+params.ρ_caco3)*P+params.α_p*params.γ*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*params.Rd_phy*P+params.α_z*params.μ_z*params.Rd_phy*Z+params.α_d*params.μ_d*params.Rd_phy*D+params.α_dd*params.μ_dd*params.Rd_phy*DD+params.μ_dom*DOM*params.Rd_dom #+ delta(z-grid.zᵃᵃᶜ[Nz])*air_sea_flux(16, 2.002e-3, 2.311e-3, 8.0)/grid.Δzᵃᵃᶜ
ALK_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, PAR::AbstractFloat, params) = params.μ_p*L_I(PAR, params)*L_NO₃(NO₃, NH₄, params)*P-2*params.ρ_caco3*params.μ_p*L_I(PAR, params)*(L_NO₃(NO₃, NH₄, params)+L_NH₄(NH₄, params))*params.Rd_phy*P

OXY_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, PAR::AbstractFloat, params) = nothing

getPAR(PAR::Field, x, y, z, t) = Oceananigans.Fields.interpolate(PAR, Center(), Center(), Center(), PAR.grid, x, y, z)
getPAR(PAR, x, y, z, t) = PAR(x, y, z, t)

#source functions with PAR func
Z_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, params) = Z_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, getPAR(params.PAR, x, y, z, t), params)
D_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, params) = D_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, getPAR(params.PAR, x, y, z, t), params)
DD_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, params) = DD_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, getPAR(params.PAR, x, y, z, t), params)
P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, params) =P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, getPAR(params.PAR, x, y, z, t), params)
NO₃_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, params) = NO₃_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, getPAR(params.PAR, x, y, z, t), params)
NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, params) = NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, getPAR(params.PAR, x, y, z, t), params)
DOM_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, params) = DOM_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, getPAR(params.PAR, x, y, z, t), params)

DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)
ALK_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = ALK_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, getPAR(params.PAR, x, y, z, t), params)

OXY_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, params) = nothing

#advective forcings
D_sinking(z, params) = params.V_d*tanh(max(-z/params.λ,0))
DD_sinking(z, params) = params.V_dd*tanh(max(-z/params.λ,0))

#tracers
tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM)
optional_tracers=(carbonates=(:DIC, :ALK), oxygen=(:OXY,))

forcing_functions=(NO₃=NO₃_forcing, NH₄=NH₄_forcing, P=P_forcing, Z=Z_forcing, D=D_forcing, DD=D_forcing, DOM=DOM_forcing, DIC=DIC_forcing, ALK=ALK_forcing)
sinking=(D=D_sinking, DD=DD_sinking)
end # module

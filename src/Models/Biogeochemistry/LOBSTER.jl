"
LOBSTER biogeochemistry model as described by Karleskind et al. 2011 and earlier publications.

Variables: NO₃, NH₄, P, Z, D, DD, DOM
Optional variables: (DIC, ALK), O₂
Sinking variables: D, DD

Forcing: PAR

References:
Karleskind, P., Lévy, M., and Memery, L. (2011), Subduction of carbon, nitrogen, and oxygen in the northeast Atlantic, J. Geophys. Res., 116, C02025, doi:10.1029/2010JC006446. 
"
module LOBSTER
#https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2004JC002588
using Oceananigans
using Oceananigans.Units: second, minute, minutes, hour, hours, day, days, year, years

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

OXY_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, PAR::AbstractFloat, params) = params.μ_p*L_I(PAR, params)*L_NO₃(NO₃, NH₄, params)*params.Rd_oxy*P-(params.Rd_oxy-params.Rd_nit)*NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, PAR, params) - params.Rd_oxy*params.μ_n*NH₄#there is a typo in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2010JC006446 so I am not sure the first term of this is correct, but this makes sense

#=
doesn't work with nonlinear grid I think (??)
function getPAR(PAR::Field, x, y, z, t)
    #assuming both are center fields, otherwise need to revert to interpolations
    #new method is about 1/3 faster
    #Oceananigans.Fields.interpolate(PAR, Center(), Center(), Center(), PAR.grid, x, y, z)
    (fi, i::Int), (fj, j::Int), (fk, k::Int) = modf.(Oceananigans.Fields.fractional_indices(x, y, z, (Center(), Center(), Center()), PAR.grid))
    return PAR[i+1, j+1, k+1]
end=#
getPAR(PAR::Field, x, y, z, t) = Oceananigans.Fields.interpolate(PAR, Center(), Center(), Center(), PAR.grid, x, y, round(z, digits=5))
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

OXY_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, params) = OXY_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, getPAR(params.PAR, x, y, z, t), params)

#advective forcings
D_sinking(z, params) = params.V_d*tanh(max(-z/params.λ,0))
DD_sinking(z, params) = params.V_dd*tanh(max(-z/params.λ,0))

#tracers
tracers=(:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM)
optional_tracers=(carbonates=(:DIC, :ALK), oxygen=(:OXY,))

forcing_functions=(NO₃=NO₃_forcing, NH₄=NH₄_forcing, P=P_forcing, Z=Z_forcing, D=D_forcing, DD=DD_forcing, DOM=DOM_forcing, DIC=DIC_forcing, ALK=ALK_forcing, OXY=OXY_forcing)
sinking=(D=D_sinking, DD=DD_sinking)

#parameters
const defaults = (
    p̃ = 0.5,  # Preference for phytoplankton
    g_z = 9.26e-6,   # Zooplankton maximal grazing rate  s⁻¹
    K_z = 1.0,    # Grazing half-saturation value    mmolm⁻³
    #Q_sol = 600,  # incoming solar radiation   Wm⁻², need to check the number
    k_r0 = 0.225,  # m⁻¹
    k_b0 = 0.0232,  # m⁻¹
    Χ_rp = 0.037,  # m⁻¹(mgChlm⁻³)⁻ᵉʳ
    Χ_bp = 0.074,  # m⁻¹(mgChlm⁻³)⁻ᵉᵇ
    e_r = 0.629, 
    e_b = 0.674, 
    #R_C2N = 6.625,  #   ratio   mmoleC/mmoleN
    #R_C2Chl = 60,  #     ratio   mmoleC/mgChl
    r_pig = 0.7,
    K_par = 33,  # Light limitation half-saturation value  Wm⁻²
    ψ = 3,  # Inhibition of nitrate uptake by ammonium
    K_no₃ = 0.7,  # Nitrate limitation half-saturation value   mmolm⁻³
    K_nh₄ = 0.001,  # Ammonium limitation half-saturation value   mmolm⁻³
    v_dd_min = 50/day, #DD min sinking speed  50/day
    v_dd_max = 200/day, #DD max sinking speed
    #detritus sinking parameters
    V_d = -3.47e-5,  #  Detritus sedimentation speed   ms⁻¹
    V_dd = -200/day,  #  Detritus sedimentation speed  -v_dd_min       50m/day=0.0005878  ms⁻¹
    λ = 1, # I think it should be deeper than the first grid point     
    μ_p = 1.21e-5, #  s⁻¹   Phytoplankton maximal growth rate   1/day

    #"Lobster_parameters"
    a_z = 0.7,  # Assimilated food fraction by zooplankton
    m_z = 2.31e-6,  # Zooplankton mortality rate  s⁻¹mmol⁻¹m³
    μ_z = 5.8e-7, # Zooplankton excretion rate  s⁻¹
    #g_z = 9.26e-6, # Zooplankton maximal grazing rate  s⁻¹
    #p̃ = 0.5, # Preference for phytoplankton
    #K_z = 1.0, # Grazing half-saturation value    mmolm⁻³
    m_p = 5.8e-7, # Phytoplankton mortality rate   s⁻¹

    μ_d = 5.78e-7, #  Detritus remineralization rate  s⁻¹
    μ_dd = 5.78e-7, 
    γ = 0.05,  #  Phytoplankton exudation rate  

    μ_n = 5.8e-7, # Nitrification rate   s⁻¹
    #f_n = 0.75, #  Ammonium/DOM redistribution ratio
    α_p = 0.75,  #NH4 fraction of P exsudation
    α_z = 0.5,  #NH4 fraction of Z excretion
    α_d = 0, #NH4 fraction of D degradation
    α_dd = 0, #NH4 fraction of  DD degradation

    Rd_phy = 6.56, # C:N ratio for P molC mol N -1
    Rd_dom = 6.56,
    Rd_chl = 1.31,
    ρ_caco3 = 0.1, # rain ratio of organic carbon to CaCO3 
    Rd_oxy = 10.75, #O:N for PP
    Rd_nit = 2.0, #O:N for Nitrification

    f_z = 0.5, #  Fraction of slow sinking mortality   0.5  
    f_d = 0.5, # Faecal pellets and P mortality fraction to DD   0.5

    μ_dom = 3.86e-7, # DOM breakdown rate    s⁻¹
)
end # module

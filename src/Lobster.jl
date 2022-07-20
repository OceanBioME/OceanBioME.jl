module Lobster

using Oceananigans
using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years
using Roots
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Architectures: device

include("parameters.jl")

@inline p(P, D, params) = params.p̃*P/(params.p̃*P+(1-params.p̃)*D)
@inline G_d(P, Z, D, params) = params.g_z*(1-p(P, D, params))*D*Z/(params.K_z+P*p(P, D, params)+(1-p(P, D, params))*D)
@inline G_p(P, Z, D, params) = params.g_z*p(P, D, params)*P*Z/(params.K_z+P*p(P, D, params)+(1-p(P, D, params))*D)

#Limiting equations
@inline L_I(PAR, params) = 1-exp(-PAR/params.K_par) # reference Levy 2001
@inline L_NO₃(NO₃, NH₄, params) = NO₃*exp(-params.ψ*NH₄)/(NO₃+params.K_no₃)
@inline L_NH₄(NH₄, params) = max(NH₄/(NH₄+params.K_nh₄),0)    # NH₄/(NH₄+K_nh₄)

function air_sea_flux(x, y, t, DIC, ALK, T::AbstractFloat, S::AbstractFloat, params) # has to be only including x,y,t without z, because this will apply to the z direction. f(x, y, t) on z-boundaries.
    #https://clima.github.io/OceananigansDocumentation/stable/model_setup/boundary_conditions/
    #https://biocycle.atmos.colostate.edu/shiny/carbonate/
    #ALK *= 1.e-6 # microequivalents to equivalents
    #DIC *= 1.e-6 # micromoles to moles
    ALK *= 1.e-3/params.ρₒ # microequivalents to equivalents  from mmol/m^-3 to mol/kg
    DIC *= 1.e-3/params.ρₒ # micromoles to moles    

    Boron = 1.179e-5*S # Total Boron mole/kg as a fraction of salinity
    
    K0 = exp(-60.2409 + 9345.17/T + 23.3585*log(T/100) + S*(0.023517 - 0.00023656*T + 0.0047036*(T/100)^2))   # mol/kg/atm 
    K1 = exp(2.18867 - 2275.036/T - 1.468591*log(T) + (-0.138681 - 9.33291/T)*sqrt(S) + 0.0726483*S - 0.00574938*S^1.5)
    K2 = exp(-0.84226 - 3741.1288/T -1.437139*log(T) + (-0.128417 - 24.41239/T)*sqrt(S) + 0.1195308*S - 0.0091284*S^1.5)
    KB = exp( (-8966.90 - 2890.51*sqrt(S) - 77.942*S + 1.726*S^1.5 - 0.0993*S^2)/T + (148.0248 + 137.194*sqrt(S) + 1.62247*S) + (-24.4344 - 25.085*sqrt(S) - 0.2474*S)*log(T) + 0.053105*sqrt(S)*T)

    H = 10^(-params.pH) # initial guess from arg list

    #potential peformance bott;eneck
    @inline CA(H) = ALK - (KB/(KB + H))*Boron
    @inline H_eq(H) = CA(H)*H^2 + K1*(CA(H)-DIC)*H + K1*K2*(CA(H)-2*DIC)
    H = find_zero(H_eq, H)
    
    #pH = -log10(H)
    CO2aq = CA(H)/(K1/H + 2*K1*K2/H^2)*1e6 # Eq 11  μmol/kg
    pCO2 = CO2aq/K0 # Eq 4 (converted from atm to ppmv) ppm 
    HCO3 = K0*K1*CO2aq/K0/H       #μmol/kg
    CO3 = K0*K1*K2*CO2aq/K0/H^2   #μmol/kg
    DIC = CO2aq + HCO3 + CO3
    R = DIC/CO3
    #[flux, pCO2, CO2aq, HCO3, CO3, DIC, R]
    flux = 7.7e-4*params.U_10^2*(pCO2-params.pCO2_air)/(365*24*3600/1000) # mmol/m^2/S

    return flux      #Wanninkhof 2014 equ.6 positive value means upward flux meaning losing Carbon 
end

function air_sea_flux(x, y, t, DIC, ALK, params)
    return air_sea_flux(x, y, t, DIC, ALK, params.T(x, y, 0, t), params.S(x, y, 0, t), params)
end

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

Z_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = Z_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)
D_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = D_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)
DD_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = DD_forcing_nonslip(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)
P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) =P_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)
NO₃_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = NO₃_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)
NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = NH₄_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)
DOM_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = DOM_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)
DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = DIC_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)
ALK_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params) = ALK_forcing(x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, params.PAR(x, y, z, t), params)

mutable struct BGCModel{T, F, B}
    tracers :: T
    forcing :: F
    boundary_conditions :: B
end

function lobster(grid, parameters, forcings=(T=nothing, S=nothing, PAR=nothing))
    if keys(forcings) != (:T, :S, :PAR)
        throw(ArgumentError("Need to provide all external forcings (T, S, PAR), either pass 'nothing' for a field to be used, or pass a function of x, y, z, t"))
    end
    #setup forcings
    slip_vel_D = zeros(0:grid.Nx+2,0:grid.Ny+2,0:grid.Nz+2)
    slip_vel_DD = zeros(0:grid.Nx+2,0:grid.Ny+2,0:grid.Nz+2)
    @simd for k=0:grid.Nz-2           #0:Nz+2   3:Nz-2
        #slip_vel_D[:,:,k].+=V_d*(tanh(max(-grid.zᵃᵃᶠ[k]/λ,0)))#*tanh(max((grid.zᵃᵃᶠ[k]+Lz)/λ,0)));
        @inbounds slip_vel_D[:,:,k].+=parameters.V_d*tanh(max(-grid.zᵃᵃᶠ[k]/parameters.λ,0))#*tanh(max((grid.zᵃᵃᶠ[k]+Lz)/λ,0));
        @inbounds slip_vel_DD[:,:,k].+=parameters.V_dd*tanh(max(-grid.zᵃᵃᶠ[k]/parameters.λ,0))#*tanh(max((grid.zᵃᵃᶠ[k]+Lz)/λ,0));
    end
    D_RHS_slip = AdvectiveForcing(WENO5(; grid), w=slip_vel_D)
    DD_RHS_slip = AdvectiveForcing(WENO5(; grid), w=slip_vel_DD)

    if isnothing(forcings.T)
        function_field_dependencies = (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM, :DIC, :ALK, :PAR)
    else 
        function_field_dependencies = (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM, :DIC, :ALK)
        function_parameters = merge(parameters, (PAR = forcings.PAR, ))
    end

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

    airseaflux_bc = FluxBoundaryCondition(air_sea_flux, field_dependencies = bcs_field_dependencies, parameters = bcs_parameters)
    DIC_bcs = FieldBoundaryConditions(top = airseaflux_bc, bottom = FluxBoundaryCondition(0))             #airseaflux_bc, #FluxBoundaryCondition(-1.3e-4)                                
    ALK_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(0), bottom = FluxBoundaryCondition(0))

    lobster = BGCModel((:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM, :DIC, :ALK), #tracers
        (Z=Z_RHS, D=(D_RHS_nonslip, D_RHS_slip), DD=(DD_RHS_nonslip, DD_RHS_slip), P=P_RHS, NO₃=NO₃_RHS, NH₄=NH₄_RHS, DOM=DOM_RHS, DIC=DIC_RHS, ALK=ALK_RHS), #forcing
        (P=P_bcs, Z=Z_bcs, D=D_bcs, DD=DD_bcs, NO₃=NO₃_bcs, NH₄=NH₄_bcs, DOM=DOM_bcs, DIC=DIC_bcs, ALK=ALK_bcs)) #boundaries

    return lobster
end

function setup(grid, parameters, tracers = (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM, :DIC, :ALK)) #could here use the tracers argument to use different models, e.g. tracers = (:N, :P, :Z) and then change this setup function to return different functions
    if tracers != (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM, :DIC, :ALK)
        throw(ArgumentError("Only tracer combination implimented is  (:NO₃, :NH₄, :P, :Z, :D, :DD, :DOM, :DIC, :ALK)"))
    else
        return lobster(grid, parameters)
    end
end

@kernel function _update_par!(par, grid, chl, t, params) 
    #from https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/JC093iC09p10749 and https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2000JC000319
    i, j = @index(Global, NTuple) 
 
    ∫chl = chl[i, j, grid.Nz] * Oceananigans.Operators.Δzᶜᶜᶜ(i, j, grid.Nz, grid) 
    beta = 0.0
    @unroll for band in [1:length(params.kw_bands);]
        beta += params.kw_bands[band]+params.χ_bands[band]*∫chl^params.e_bands[band]
    end
    @inbounds par[i, j, grid.Nz] = params.surface_PAR(mod(t, 364days))*exp(beta*grid.zᵃᵃᶜ[grid.Nz])
 
    #@unroll for k in grid.Nz-1 : -1 : 1 
    for k in grid.Nz-1 : -1 : 1 
        @inbounds ∫chl += chl[i, j, k+1] * Oceananigans.Operators.Δzᶜᶜᶜ(i, j, k+1, grid) 
        beta = 0.0
        #@unroll for band in [1:length(params.kw_bands)]
        for band in [1:length(params.kw_bands);]
            beta += params.kw_bands[band]+params.χ_bands[band]*∫chl^params.e_bands[band]
        end
        @inbounds par[i, j, k] = params.surface_PAR(mod(t, 364days))*exp(beta*grid.zᵃᵃᶜ[k])
    end 
end 
function update_par!(sim, params)
    par_calculation = Oceananigans.Utils.launch!(sim.model.architecture, sim.model.grid, :xy, _update_par!,
                                   sim.model.tracers.PAR, sim.model.grid, sim.model.tracers.P .* params.Rd_phy, time(sim), params,
                                   dependencies = Event(device(sim.model.architecture)))

    wait(device(sim.model.architecture), par_calculation)
end
end # module

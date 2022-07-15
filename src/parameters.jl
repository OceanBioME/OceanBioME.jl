using Oceananigans.Units: second,minute, minutes, hour, hours, day, days, year, years

const default = (
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
    #r_pig = 0.7,
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
    ρ_caco3 = 0.1, # rain ratio of organic carbon to CaCO3 

    f_z = 0.5, #  Fraction of slow sinking mortality   0.5  
    f_d = 0.5, # Faecal pellets and P mortality fraction to DD   0.5

    pH=8.0, # initial pH value guess for air-sea flux calculation
    U_10 = 10.0,  # 10m/s for CO2 flux calculation only
    pCO2_air = 413.3, #  in Jan 2020   https://www.co2.earth/  for CO2 flux calculation only
    μ_dom = 3.86e-7, # DOM breakdown rate    s⁻¹

    ρₒ = 1026 # kg m⁻³, average density at the surface of the world ocean
)
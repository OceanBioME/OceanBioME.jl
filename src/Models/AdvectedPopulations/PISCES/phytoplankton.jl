include("base_production.jl")
include("nutrient_limitation.jl")

@kwdef struct Phytoplankton{GR, FT}
                        growth_rate :: GR
                        
                   exudated_fracton :: FT = 0.05

              blue_light_absorption :: FT
             green_light_absorption :: FT
               red_light_absorption :: FT

          mortality_half_saturation :: FT = 0.2
              linear_mortality_rate :: FT = 0.01

           base_quadratic_mortality :: FT = 0.01
        maximum_quadratic_mortality :: FT # zero for nanophytoplankton

          minimum_chlorophyll_ratio :: FT = 0.0033
          maximum_chlorophyll_ratio :: FT

                 maximum_iron_ratio :: FT = 40.0

           silicate_half_saturation :: FT = 2.0
  enhanced_silicate_half_saturation :: FT = 20.9
             optimal_silicate_ratio :: FT = 0.159
end

@inline phytoplankton_concentration(::Val{P}, P, D) = P
@inline phytoplankton_concentration(::Val{D}, P, D) = D

@inline phytoplankton_grazing(::Val{P}, args...) = nanophytoplankton_grazing(args...)
@inline phytoplankton_grazing(::Val{D}, args...) = diatom_grazing(args...)

@inline function (phyto::Phytoplankton)(bgc, val_name::Union{Val{:P}, Val{:D}}, 
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi, # we should get rid of DSi and the rest of the Si since it doesn't do anything...
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, 
                                        zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)
    # production
    δ  = phyto.exudated_fracton
    β₁ = phyto.blue_light_absorption
    β₂ = phyto.green_light_absorption
    β₃ = phyto.red_light_absorption

    PARᵢ = β₁ * PAR₁ + β₂ * PAR₂ + β₃ * PAR₃

    I    = phytoplankton_concentration(val_name, P, D)
    IChl = phytoplankton_concentration(val_name, PChl, DChl)
    IFe  = phytoplankton_concentration(val_name, PFe, DFe)

    L, = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μ = phyto.growth_rate(bgc, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PARᵢ, L)

    production = (1 - δ) * μ * I 

    # mortality
    K = phyto.mortality_half_saturation
    m = phyto.linear_mortality_rate

    backgroound_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear

    linear_mortality = m * I / (I + K) * I

    w₀ = phyto.base_quadratic_mortality
    w₁ = phyto.maximum_quadratic_mortality

    w = w₀ + w₁ * (1 - L)
    
    shear = ifelse(z < zₘₓₗ, backgroound_shear, mixed_layer_shear)

    quadratic_mortality = shear * w * I^2

    # grazing
    gZ = phytoplankton_grazing(val_name, bgc.zooplankton, P, D, Z, POC)
    gM = phytoplankton_grazing(val_name, bgc.mesozooplankton, P, D, Z, POC)

    grazing = gZ * Z + gM * M

    return production - linear_mortality - quadratic_mortality - grazing
end

@inline function (phyto::Phytoplankton)(bgc, val_name::Union{Val{:PChl}, Val{:DChl}}, 
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi, # we should get rid of DSi and the rest of the Si since it doesn't do anything...
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, 
                                        zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)

    I    = phytoplankton_concentration(val_name, P, D)
    IChl = phytoplankton_concentration(val_name, PChl, DChl)
    IFe  = phytoplankton_concentration(val_name, PFe, DFe)

    # production
    δ  = phyto.exudated_fracton
    β₁ = phyto.blue_light_absorption
    β₂ = phyto.green_light_absorption
    β₃ = phyto.red_light_absorption

    θ₀ = phyto.minimum_chlorophyll_ratio
    θ₁ = phyto.maximum_chlorophyll_ratio

    PARᵢ = β₁ * PAR₁ + β₂ * PAR₂ + β₃ * PAR₃
    L = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μ, ρ = production_and_energy_assimilation_absorption_ratio(phyto.growth_rate, bgc, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PARᵢ, L)

    production = (1 - δ) * (12 * θ₀ + (θ₁ - θ₀) * ρ) * μ * I 

    # mortality
    K = phyto.mortality_half_saturation
    m = phyto.linear_mortality_rate

    backgroound_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear

    linear_mortality = m * I / (I + K) * IChl

    w₀ = phyto.base_quadratic_mortality
    w₁ = phyto.maximum_quadratic_mortality

    w = w₀ + w₁ * (1 - L)
    
    shear = ifelse(z < zₘₓₗ, backgroound_shear, mixed_layer_shear)

    quadratic_mortality = shear * w * I * IChl

    # grazing
    θChl = IChl / I

    gZ = phytoplankton_grazing(val_name, bgc.zooplankton, P, D, Z, POC)
    gM = phytoplankton_grazing(val_name, bgc.mesozooplankton, P, D, Z, POC)

    grazing = (gZ * Z + gM * M) * θChl

    return production - linear_mortality - quadratic_mortality - grazing
end

@inline function (phyto::Phytoplankton)(bgc, val_name::Union{Val{:PFe}, Val{:DFe}}, 
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi, # we should get rid of DSi and the rest of the Si since it doesn't do anything...
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, 
                                        zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)

    I    = phytoplankton_concentration(val_name, P, D)
    IChl = phytoplankton_concentration(val_name, PChl, DChl)
    IFe  = phytoplankton_concentration(val_name, PFe, DFe)

    # production
    δ  = phyto.exudated_fracton

    θFeₘ = phyto.maximum_iron_ratio

    θFe = IFe / I

    L, LFe = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μᵢ = base_production_rate(phyto.growth_rate, bgc, T)

    L₁ = iron_uptake_limitation(phyto.growth_rate, I, Fe) # assuming bFe = Fe

    L₂ = (4 - 2 * LFe) / (LFe + 1) # Formulation in paper does not vary between 1 and 4 as claimed, this does

    production = (1 - δ) * θFeₘ * L₁ * L₂ * (1 - θFe / θFeₘ) / (1.05 - θFe / θFeₘ) * μᵢ * I 

    # mortality
    K = phyto.mortality_half_saturation
    m = phyto.linear_mortality_rate

    backgroound_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear

    linear_mortality = m * I / (I + K) * IFe

    w₀ = phyto.base_quadratic_mortality
    w₁ = phyto.maximum_quadratic_mortality

    w = w₀ + w₁ * (1 - L)
    
    shear = ifelse(z < zₘₓₗ, backgroound_shear, mixed_layer_shear)

    quadratic_mortality = shear * w * I * IFe

    # grazing
    gZ = phytoplankton_grazing(val_name, bgc.zooplankton, P, D, Z, POC)
    gM = phytoplankton_grazing(val_name, bgc.mesozooplankton, P, D, Z, POC)

    grazing = (gZ * Z + gM * M) * θFe

    return production - linear_mortality - quadratic_mortality - grazing
end

@inline function (phyto::Phytoplankton)(bgc, ::Val{:DSi}, 
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi, # we should get rid of DSi and the rest of the Si since it doesn't do anything...
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, 
                                        zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)

    # production
    δ  = phyto.exudated_fracton
    β₁ = phyto.blue_light_absorption
    β₂ = phyto.green_light_absorption
    β₃ = phyto.red_light_absorption

    K₁ = phyto.silicate_half_saturation
    K₂ = phyto.enhanced_silicate_half_saturation
    θ₀ = phyto.optimal_silicate_ratio

    PARᵢ = β₁ * PAR₁ + β₂ * PAR₂ + β₃ * PAR₃

    φ = bgc.latitude(y)

    L, LFe, LPO₄ = phyto.nutrient_limitation(bgc, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μ = phyto.growth_rate(bgc, D, DChl, T, zₘₓₗ, zₑᵤ, κ, PARᵢ, L)
    μᵢ = base_production_rate(phyto.growth_rate, bgc, T)

    L₁ = Si / (Si + K₁)

    # enhanced silication in southern ocean
    L₂ = ifelse(φ < 0, Si^3 / (Si^3 + K₂^3), 0)

    F₁ = min(μ / (μᵢ * L), LFe, LPO₄, LN)

    F₂ = min(1, 2.2 * max(0, L₁ - 0.5))

    θ₁ = θ₀ * L₁ * min(5.4, (4.4 * exp(-4.23 * F₁) * F₂ + 1) * (1 + 2 * L₂))

    production = (1 - δ) * θ₁ * μ * D

    # mortality
    K = phyto.mortality_half_saturation
    m = phyto.linear_mortality_rate

    backgroound_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear

    linear_mortality = m * D / (D + K) * IFe

    w₀ = phyto.base_quadratic_mortality
    w₁ = phyto.maximum_quadratic_mortality

    w = w₀ + w₁ * (1 - L)
    
    shear = ifelse(z < zₘₓₗ, backgroound_shear, mixed_layer_shear)

    quadratic_mortality = shear * w * D * IFe

    # grazing
    gZ = diatom_grazing(bgc.zooplankton, P, D, Z, POC)
    gM = diatom_grazing(bgc.mesozooplankton, P, D, Z, POC)

    grazing = (gZ * Z + gM * M) * θFe

    return production - linear_mortality - quadratic_mortality - grazing
end

include("base_production.jl")
include("nutrient_limitation.jl")

@kwdef struct Phytoplankton{GR, NL, FT}
                        growth_rate :: GR
                nutrient_limitation :: NL
                        
                   exudated_fracton :: FT = 0.05

              blue_light_absorption :: FT
             green_light_absorption :: FT
               red_light_absorption :: FT

          mortality_half_saturation :: FT = 0.2
              linear_mortality_rate :: FT = 0.01 / day

           base_quadratic_mortality :: FT = 0.01 / day
        maximum_quadratic_mortality :: FT # zero for nanophytoplankton

          minimum_chlorophyll_ratio :: FT = 0.0033
          maximum_chlorophyll_ratio :: FT

                 maximum_iron_ratio :: FT = 40.0e-3

           silicate_half_saturation :: FT = 2.0
  enhanced_silicate_half_saturation :: FT = 20.9
             optimal_silicate_ratio :: FT = 0.159
end

@inline phytoplankton_concentration(::NANO_PHYTO, P, D) = P
@inline phytoplankton_concentration(::DIATOMS, P, D) = D

@inline phytoplankton_grazing(::NANO_PHYTO, args...) = nanophytoplankton_grazing(args...)
@inline phytoplankton_grazing(::DIATOMS, args...) = diatom_grazing(args...)

@inline function (phyto::Phytoplankton)(val_name::Union{Val{:P}, Val{:D}}, bgc,
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi, # we should get rid of DSi and the rest of the Si since it doesn't do anything...
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, S,
                                        zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)
    # production
    δ  = phyto.exudated_fracton

    I    = phytoplankton_concentration(val_name, P, D)
    IChl = phytoplankton_concentration(val_name, PChl, DChl)
    IFe  = phytoplankton_concentration(val_name, PFe, DFe)

    μI, L = total_production(phyto, bgc, y, t, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    production = (1 - δ) * μI 

    # mortality
    linear_mortality, quadratic_mortality = mortality(phyto, bgc, z, I, zₘₓₗ, L)

    # grazing
    gZ = phytoplankton_grazing(val_name, bgc.microzooplankton, P, D, Z, POC, T)
    gM = phytoplankton_grazing(val_name, bgc.mesozooplankton, P, D, Z, POC, T)

    grazing = gZ * Z + gM * M   

    return production - linear_mortality - quadratic_mortality - grazing
end

@inline function (phyto::Phytoplankton)(val_name::Union{Val{:PChl}, Val{:DChl}}, bgc,
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi, # we should get rid of DSi and the rest of the Si since it doesn't do anything...
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, S,
                                        zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    I    = phytoplankton_concentration(val_name, P, D)
    IChl = phytoplankton_concentration(val_name, PChl, DChl)
    IFe  = phytoplankton_concentration(val_name, PFe, DFe)

    # production
    δ  = phyto.exudated_fracton

    θ₀ = phyto.minimum_chlorophyll_ratio
    θ₁ = phyto.maximum_chlorophyll_ratio

    L, = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μ, ρ = production_and_energy_assimilation_absorption_ratio(phyto.growth_rate, phyto, bgc, y, t, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR, PAR₁, PAR₂, PAR₃, L)

    production = (1 - δ) * (12 * θ₀ + (θ₁ - θ₀) * ρ) * μ * I 

    # mortality
    θChl = IChl / (12 * I + eps(0.0))

    linear_mortality, quadratic_mortality = mortality(phyto, bgc, z, I, zₘₓₗ, L)

    linear_mortality *= θChl
    quadratic_mortality *= θChl
    
    # grazing

    gZ = phytoplankton_grazing(val_name, bgc.microzooplankton, P, D, Z, POC, T)
    gM = phytoplankton_grazing(val_name, bgc.mesozooplankton, P, D, Z, POC, T)

    grazing = (gZ * Z + gM * M) * θChl

    return production - linear_mortality - quadratic_mortality - grazing
end

@inline function (phyto::Phytoplankton)(val_name::Union{Val{:PFe}, Val{:DFe}}, bgc,
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi, # we should get rid of DSi and the rest of the Si since it doesn't do anything...
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, S,
                                        zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    I    = phytoplankton_concentration(val_name, P, D)
    IChl = phytoplankton_concentration(val_name, PChl, DChl)
    IFe  = phytoplankton_concentration(val_name, PFe, DFe)

    # production
    production, L = iron_uptake(phyto, bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T)

    # mortality
    linear_mortality, quadratic_mortality = mortality(phyto, bgc, z, I, zₘₓₗ, L)

    linear_mortality *= IFe / (I + eps(0.0))
    quadratic_mortality *= IFe / (I + eps(0.0))
    
    # grazing
    gZ = phytoplankton_grazing(val_name, bgc.microzooplankton, P, D, Z, POC, T)
    gM = phytoplankton_grazing(val_name, bgc.mesozooplankton, P, D, Z, POC, T)

    grazing = (gZ * Z + gM * M) * IFe / (I + eps(0.0))

    return production - linear_mortality - quadratic_mortality - grazing
end

@inline function iron_uptake(phyto::Phytoplankton, bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T)
    δ  = phyto.exudated_fracton
    θFeₘ = phyto.maximum_iron_ratio

    θFe = IFe / (I + eps(0.0))

    L, LFe = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μᵢ = base_production_rate(phyto.growth_rate, T)

    L₁ = iron_uptake_limitation(phyto.nutrient_limitation, I, Fe) # assuming bFe = Fe

    L₂ = (4 - 2 * LFe) / (LFe + 1) # Formulation in paper does not vary between 1 and 4 as claimed, this does

    return (1 - δ) * θFeₘ * L₁ * L₂ * (1 - θFe / θFeₘ) / (1.05 - θFe / θFeₘ) * μᵢ * I, L
end

@inline function (phyto::Phytoplankton)(::Val{:DSi}, bgc,
                                        x, y, z, t,
                                        P, D, Z, M, 
                                        PChl, DChl, PFe, DFe, DSi,
                                        DOC, POC, GOC, 
                                        SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, 
                                        CaCO₃, DIC, Alk, 
                                        O₂, T, S,
                                        zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    # production
    production, L = silicate_uptake(phyto, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    # mortality
    linear_mortality, quadratic_mortality = mortality(phyto, bgc, z, D, zₘₓₗ, L)

    linear_mortality *= DSi / (D + eps(0.0))
    quadratic_mortality *= DSi / (D + eps(0.0))

    # grazing
    gZ = diatom_grazing(bgc.microzooplankton, P, D, Z, POC, T)
    gM = diatom_grazing(bgc.mesozooplankton, P, D, Z, POC, T)

    grazing = (gZ * Z + gM * M) * DSi / (D + eps(0.0))

    return production - linear_mortality - quadratic_mortality - grazing
end

@inline function silicate_uptake(phyto::Phytoplankton, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    δ  = phyto.exudated_fracton

    K₁ = phyto.silicate_half_saturation
    K₂ = phyto.enhanced_silicate_half_saturation
    θ₀ = phyto.optimal_silicate_ratio
    
    L, LFe, LPO₄, LN = phyto.nutrient_limitation(bgc, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μ = phyto.growth_rate(phyto, bgc, y, t, D, DChl, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃, L)

    μᵢ = base_production_rate(phyto.growth_rate, T)

    L₁ = Si / (Si + K₁ + eps(0.0))

    # enhanced silication in southern ocean
    φ = bgc.latitude(y)
    
    L₂ = ifelse(φ < 0, Si^3 / (Si^3 + K₂^3), 0)

    F₁ = min(μ / (μᵢ * L + eps(0.0)), LFe, LPO₄, LN)

    F₂ = min(1, 2.2 * max(0, L₁ - 0.5))

    θ₁ = θ₀ * L₁ * min(5.4, (4.4 * exp(-4.23 * F₁) * F₂ + 1) * (1 + 2 * L₂))

    return (1 - δ) * θ₁ * μ * D, L
end

@inline function dissolved_exudate(phyto::Phytoplankton, bgc, y, t, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    δ = phyto.exudated_fracton

    μI, = total_production(phyto, bgc, y, t, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    return δ * μI
end

@inline function mortality(phyto::Phytoplankton, bgc, z, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′, zₘₓₗ)
    L, = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    return mortality(phyto, bgc, z, I, zₘₓₗ, L)
end

@inline function mortality(phyto::Phytoplankton, bgc, z, I, zₘₓₗ, L)
    K = phyto.mortality_half_saturation
    m = phyto.linear_mortality_rate

    background_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear

    linear_mortality = m * I / (I + K) * I

    w₀ = phyto.base_quadratic_mortality
    w₁ = phyto.maximum_quadratic_mortality

    # this enhanced mortality is parameterised as per NEMO but differs from Aumount 2015 (w₁ * (1 - L))
    w = w₀ + w₁ * 0.25 * (1 - L^2) / (0.25 + L^2)
    
    shear = ifelse(z < zₘₓₗ, background_shear, mixed_layer_shear)

    quadratic_mortality = shear * w * I^2

    return linear_mortality, quadratic_mortality
end

@inline function nitrate_uptake(phyto::Phytoplankton, bgc, y, t, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    L, _, _, LN, L_NO₃ = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μ = phyto.growth_rate(phyto, bgc, y, t, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃, L) * I

    return μ * L_NO₃ / (LN + eps(0.0))
end

@inline function ammonia_uptake(phyto::Phytoplankton, bgc, y, t, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    L, _, _, LN, _, L_NH₄ = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    μ = phyto.growth_rate(phyto, bgc, y, t, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃, L) * I

    return μ * L_NH₄ / (LN + eps(0.0))
end

@inline function total_production(phyto::Phytoplankton, bgc, y, t, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, T, Si′, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)
    L, = phyto.nutrient_limitation(bgc, I, IChl, IFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    return phyto.growth_rate(phyto, bgc, y, t, I, IChl, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃, L) * I, L
end

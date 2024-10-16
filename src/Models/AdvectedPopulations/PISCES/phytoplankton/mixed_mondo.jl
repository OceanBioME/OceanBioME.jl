using Oceananigans.Grids: znode, Center

"""
    MixedMondo

Holds the parameters for the PISCES mixed mondo phytoplankton 
parameterisation where nutrient limitation is modelled using the
mondo approach for nitrate (NO₃), ammonia (NH₄), phosphate (PO₄),
and silicate (Si), but the quota approach is used for iron (Fe) 
and light (PAR).

Therefore each class has a carbon compartement (generically `I`),
chlorophyll (`IChl`), and iron (`IFe`), and may also have silicate
(`ISi`) if the `nutrient_limitation` specifies that the growth is
silicate limited, despite the fact that the silicate still limits 
the growth in a mondo fashion.

The `growth_rate` may be different parameterisations, currently 
either `NutrientLimitedProduction` or 
`GrowthRespirationLimitedProduction`, which represent the typical
and `newprod` versions of PISCES.
"""
struct MixedMondo{GR, NL, FT}
                        growth_rate :: GR
                nutrient_limitation :: NL
                        
                   exudated_fracton :: FT # 

              blue_light_absorption :: FT #
             green_light_absorption :: FT #
               red_light_absorption :: FT #

          mortality_half_saturation :: FT # mmol C / m³
              linear_mortality_rate :: FT # 1 / s

           base_quadratic_mortality :: FT # 1 / s / (mmol C / m³)
        maximum_quadratic_mortality :: FT # 1 / s / (mmol C / m³) - zero for nanophytoplankton

          minimum_chlorophyll_ratio :: FT # mg Chl / mg C
          maximum_chlorophyll_ratio :: FT # mg Chl / mg C

                 maximum_iron_ratio :: FT # μmol Fe / mmol C

           silicate_half_saturation :: FT # mmol Si / m³
  enhanced_silicate_half_saturation :: FT # mmol Si / m³
             optimal_silicate_ratio :: FT # mmol Si / mmol C

    half_saturation_for_iron_uptake :: FT # μmol Fe / m³

      threshold_for_size_dependency :: FT # mmol C / m³
                         size_ratio :: FT # 

    function MixedMondo(FT = Float64;
                        growth_rate::GR,
                        nutrient_limitation::NL,
                        exudated_fracton = 0.05, # 

                        blue_light_absorption, #
                        green_light_absorption, #
                        red_light_absorption, #
          
                        mortality_half_saturation = 0.2, # mmol C / m³
                        linear_mortality_rate = 0.01 / day, # 1 / s
          
                        base_quadratic_mortality = 0.01 / day, # 1 / s / (mmol C / m³)
                        maximum_quadratic_mortality, # 1 / s / (mmol C / m³) - zero for nanophytoplankton
          
                        minimum_chlorophyll_ratio = 0.0033, # mg Chl / mg C
                        maximum_chlorophyll_ratio, # mg Chl / mg C
          
                        maximum_iron_ratio = 0.06, # μmol Fe / mmol C
          
                        silicate_half_saturation = 2.0, # mmol Si / m³
                        enhanced_silicate_half_saturation = 20.9, # mmol Si / m³
                        optimal_silicate_ratio = 0.159, # mmol Si / mmol C
          
                        half_saturation_for_iron_uptake,  # μmol Fe / m³
          
                        threshold_for_size_dependency = 1.0,  # mmol C / m³
                        size_ratio = 3.0) where {GR, NL} #

        return new{GR, NL, FT}(growth_rate, nutrient_limitation, 
                               exudated_fracton, 
                               blue_light_absorption, green_light_absorption, red_light_absorption,
                               mortality_half_saturation, linear_mortality_rate,
                               base_quadratic_mortality, maximum_quadratic_mortality,
                               minimum_chlorophyll_ratio, maximum_chlorophyll_ratio,
                               maximum_iron_ratio, 
                               silicate_half_saturation, enhanced_silicate_half_saturation, optimal_silicate_ratio,
                               half_saturation_for_iron_uptake, 
                               threshold_for_size_dependency, size_ratio)
    end
end

required_biogeochemical_tracers(phyto::MixedMondo, base) =
    (base, Symbol(base, :Chl), Symbol(base, :Fe), 
     ifelse(phyto.nutrient_limitation.silicate_limited, (Symbol(base, :Si), ), ())...)

#####
##### Production/mortality functions
#####

@inline function carbon_growth(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    # production
    δ  = phyto.exudated_fracton

    μI = total_production(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return (1 - δ) * μI
end

@inline function chlorophyll_growth(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    I, IChl, IFe = phytoplankton_concentrations(val_name, i, j, k, fields)

    # production
    δ  = phyto.exudated_fracton

    θ₀ = phyto.minimum_chlorophyll_ratio
    θ₁ = phyto.maximum_chlorophyll_ratio

    μ, ρ = production_and_energy_assimilation_absorption_ratio(phyto.growth_rate, val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields)

    return (1 - δ) * 12 * (θ₀ + (θ₁ - θ₀) * ρ) * μ * I 
end

# production (we already account for the (1 - δ) term because it just goes straight back into Fe)
@inline iron_growth(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) = 
    iron_uptake(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline silicate_growth(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) = 
    silicate_uptake(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

#####
##### Underlying parameterisations
#####

@inline function mortality(phyto::MixedMondo, bgc, z, I, zₘₓₗ, L)
    K = phyto.mortality_half_saturation
    m = phyto.linear_mortality_rate

    background_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear

    linear_mortality = m * I / (I + K) * I

    w₀ = phyto.base_quadratic_mortality
    w₁ = phyto.maximum_quadratic_mortality

    w = w₀ + w₁ * 0.25 * (1 - L^2) / (0.25 + L^2)
    
    shear = ifelse(z < zₘₓₗ, background_shear, mixed_layer_shear)

    quadratic_mortality = shear * w * I^2

    return linear_mortality, quadratic_mortality
end

@inline function mortality(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    I, = phytoplankton_concentrations(val_name, i, j, k, fields)
    zₘₓₗ = @inbounds auxiliary_fields.zₘₓₗ[i, j, k]

    z = znode(i, j, k, grid, Center(), Center(), Center())

    L, = phyto.nutrient_limitation(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields)

    return mortality(phyto, bgc, z, I, zₘₓₗ, L)
end

@inline function total_production(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    I, = phytoplankton_concentrations(val_name, i, j, k, fields)

    L, = phyto.nutrient_limitation(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields)

    return phyto.growth_rate(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields, L) * I
end

@inline function iron_uptake(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    δ  = phyto.exudated_fracton
    θFeₘ = phyto.maximum_iron_ratio

    T  = @inbounds  fields.T[i, j, k]
    Fe = @inbounds fields.Fe[i, j, k]

    I, IChl, IFe = phytoplankton_concentrations(val_name, i, j, k, fields)

    θFe = IFe / (I + eps(0.0)) # μmol Fe / mmol C

    L, LFe = phyto.nutrient_limitation(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields)

    μᵢ = base_production_rate(phyto.growth_rate, T)

    L₁ = iron_uptake_limitation(phyto, I, Fe) # assuming bFe = Fe

    L₂ = 4 - 4.5 * LFe / (LFe + 1) # typo in Aumount 2015

    return (1 - δ) * θFeₘ * L₁ * L₂ * max(0, (1 - θFe / θFeₘ) / (1.05 - θFe / θFeₘ)) * μᵢ * I
end

@inline function iron_uptake_limitation(phyto, I, Fe)
    k = phyto.half_saturation_for_iron_uptake

    K = k * size_factor(phyto, I)

    return Fe / (Fe + K + eps(0.0))
end

@inline function size_factor(phyto, I)
    Iₘ  = phyto.threshold_for_size_dependency
    S   = phyto.size_ratio

    I₁ = min(I, Iₘ)
    I₂ = max(0, I - Iₘ)

    return (I₁ + S * I₂) / (I₁ + I₂ + eps(0.0))
end

@inline function silicate_uptake(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    δ  = phyto.exudated_fracton
    K₁ = phyto.silicate_half_saturation
    K₂ = phyto.enhanced_silicate_half_saturation
    θ₀ = phyto.optimal_silicate_ratio

    I, IChl, IFe = phytoplankton_concentrations(val_name, i, j, k, fields)

    T =  @inbounds fields.T[i, j, k]
    Si = @inbounds fields.Si[i, j, k]

    L, LFe, LPO₄, LN = phyto.nutrient_limitation(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields)

    μ = phyto.growth_rate(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields, L)

    μᵢ = base_production_rate(phyto.growth_rate, T)

    L₁ = Si / (Si + K₁ + eps(0.0))

    # enhanced silication in southern ocean
    φ = bgc.latitude(i, j, k, grid)
    
    L₂ = ifelse(φ < 0, Si^3 / (Si^3 + K₂^3), 0)

    F₁ = min(μ / (μᵢ * L + eps(0.0)), LFe, LPO₄, LN)

    F₂ = min(1, 2.2 * max(0, L₁ - 0.5))

    θ₁ = θ₀ * L₁ * min(5.4, (4.4 * exp(-4.23 * F₁) * F₂ + 1) * (1 + 2 * L₂))

    return (1 - δ) * θ₁ * μ * I
end

@inline function uptake(phyto::MixedMondo, val_name, ::Val{:NO₃}, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    _, _, _, LN, LNO₃ = phyto.nutrient_limitation(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields)

    μI = total_production(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return μI * LNO₃ / (LN + eps(0.0))
end

@inline function uptake(phyto::MixedMondo, val_name, ::Val{:NH₄}, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    _, _, _, LN, _, LNH₄ = phyto.nutrient_limitation(val_name, i, j, k, grid, bgc, phyto, clock, fields, auxiliary_fields)

    μI = total_production(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return μI * LNH₄ / (LN + eps(0.0))
end

@inline uptake(phyto::MixedMondo, val_name, ::Val{:Fe}, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    iron_uptake(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    
@inline function dissolved_exudate(phyto::MixedMondo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    δ  = phyto.exudated_fracton

    μI = total_production(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return δ * μI
end

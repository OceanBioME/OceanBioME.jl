@inline function small_mortality(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    P_linear_mortality, P_quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    R = rain_ratio(phyto, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    D_linear_mortality, = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return (1 - R / 2) * (P_linear_mortality + P_quadratic_mortality) + D_linear_mortality / 2
end

@inline function large_mortality(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    P_linear_mortality, P_quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    R = rain_ratio(phyto, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    D_linear_mortality, D_quadratic_mortality = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return R / 2 * (P_linear_mortality + P_quadratic_mortality) + D_linear_mortality / 2 + D_quadratic_mortality
end

@inline function small_mortality_iron(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    P_linear_mortality, P_quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    R = rain_ratio(phyto, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    D_linear_mortality, = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    P   = @inbounds   fields.P[i, j, k]
    PFe = @inbounds fields.PFe[i, j, k]
    D   = @inbounds   fields.D[i, j, k]
    DFe = @inbounds fields.DFe[i, j, k]

    θP = PFe / (P + eps(0.0))
    θD = DFe / (D + eps(0.0))

    return (1 - R / 2) * (P_linear_mortality + P_quadratic_mortality) * θP + D_linear_mortality * θD / 2
end

@inline function large_mortality_iron(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    P_linear_mortality, P_quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    R = rain_ratio(phyto, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    D_linear_mortality, D_quadratic_mortality = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    P   = @inbounds   fields.P[i, j, k]
    PFe = @inbounds fields.PFe[i, j, k]
    D   = @inbounds   fields.D[i, j, k]
    DFe = @inbounds fields.DFe[i, j, k]

    θP = PFe / (P + eps(0.0))
    θD = DFe / (D + eps(0.0))

    return R / 2 * (P_linear_mortality + P_quadratic_mortality) * θP + (D_linear_mortality / 2 + D_quadratic_mortality) * θD
end

@inline function coccolithophore_nutrient_limitation(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    _, _, L_PO₄, LN = phyto.nano.nutrient_limitation(Val(:P), i, j, k, grid, bgc, phyto.nano, clock, fields, auxiliary_fields)

    Fe = @inbounds fields.Fe[i, j, k]

    L_Fe = Fe / (Fe + 0.05)

    # from NEMO we should replace LPO₄ with : zlim2  = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + concnnh4 )
    # but that has to be a typo

    return min(LN, L_Fe, L_PO₄)
end

@inline coccolithophore_phytoplankton_factor(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    @inbounds max(one(grid), fields.P[i, j, k] / 2)

@inline function particulate_silicate_production(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    D   = @inbounds   fields.D[i, j, k]
    DSi = @inbounds fields.DSi[i, j, k]
    
    θ = DSi / (D + eps(0.0))

    D_linear_mortality, D_quadratic_mortality = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    total_grazing = grazing(bgc.zooplankton, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return (total_grazing + D_linear_mortality + D_quadratic_mortality) * θ
end

@inline function calcite_production(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    R = rain_ratio(phyto, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    linear_mortality, quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    total_grazing_loss = calcite_loss(bgc.zooplankton, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return R * (total_grazing_loss + (linear_mortality + quadratic_mortality) / 2)
end

@inline function rain_ratio(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    r = phyto.base_rain_ratio

    T = @inbounds fields.T[i, j, k]
    PAR = @inbounds auxiliary_fields.PAR[i, j, k]
    zₘₓₗ = @inbounds auxiliary_fields.zₘₓₗ[i, j, k]

    L_CaCO₃ = coccolithophore_nutrient_limitation(phyto, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    phytoplankton_concentration_factor = 
        coccolithophore_phytoplankton_factor(phyto, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    low_light_factor = max(0, PAR - 1) / (4 + PAR)
    high_light_factor = 30 / (30 + PAR)

    # modified from origional as it goes negative and does not achieve goal otherwise
    low_temperature_factor = max(0, T / (T + 0.1)) 
    high_temperature_factor = 1 + exp(-(T - 10)^2 / 25)

    depth_factor = min(1, -50/zₘₓₗ)

    return (r * L_CaCO₃
            * phytoplankton_concentration_factor 
            * low_light_factor 
            * high_light_factor 
            * low_temperature_factor 
            * high_temperature_factor 
            * depth_factor)
end
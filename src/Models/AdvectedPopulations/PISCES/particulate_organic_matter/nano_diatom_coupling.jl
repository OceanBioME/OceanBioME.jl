@inline function small_mortality(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields)
    P_linear_mortality, P_quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields)

    R = rain_ratio(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields)

    D_linear_mortality, = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields)

    return (1 - R / 2) * (P_linear_mortality + P_quadratic_mortality) + D_linear_mortality / 2
end

@inline function large_mortality(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields)
    P_linear_mortality, P_quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields)

    R = rain_ratio(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields)

    D_linear_mortality, D_quadratic_mortality = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields)

    return R / 2 * (P_linear_mortality + P_quadratic_mortality) + D_linear_mortality / 2 + D_quadratic_mortality
end

@inline function small_mortality_iron(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields)
    P_linear_mortality, P_quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields)

    R = rain_ratio(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields)

    D_linear_mortality, = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields)

    P   = @inbounds   fields.P[i, j, k]
    PFe = @inbounds fields.PFe[i, j, k]
    D   = @inbounds   fields.D[i, j, k]
    DFe = @inbounds fields.DFe[i, j, k]

    θP = PFe / (P + eps(0.0))
    θD = DFe / (D + eps(0.0))

    return (1 - R / 2) * (P_linear_mortality + P_quadratic_mortality) * θP + D_linear_mortality * θD / 2
end

@inline function large_mortality_iron(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields)
    P_linear_mortality, P_quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields)

    R = rain_ratio(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields)

    D_linear_mortality, D_quadratic_mortality = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields)

    P   = @inbounds   fields.P[i, j, k]
    PFe = @inbounds fields.PFe[i, j, k]
    D   = @inbounds   fields.D[i, j, k]
    DFe = @inbounds fields.DFe[i, j, k]

    θP = PFe / (P + eps(0.0))
    θD = DFe / (D + eps(0.0))

    return R / 2 * (P_linear_mortality + P_quadratic_mortality) * θP + (D_linear_mortality / 2 + D_quadratic_mortality) * θD
end

@inline function coccolithophore_nutrient_limitation(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields)
    _, _, L_PO₄, LN = phyto.nano(Val(:P), phyto.nano, i, j, k, grid, bgc, clock, fields)

    L_Fe = Fe / (Fe + 0.05)

    # from NEMO we should replace LPO₄ with : zlim2  = trb(ji,jj,jk,jppo4) / ( trb(ji,jj,jk,jppo4) + concnnh4 )
    # but that has to be a typo

    return min(LN, L_Fe, L_PO₄)
end

@inline coccolithophore_phytoplankton_factor(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields) =
    @inbounds max(one(grid), fields.P[i, j, k] / 2)

@inline function particulate_silicate_production(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields)
    D   = @inbounds   fields.D[i, j, k]
    DSi = @inbounds fields.DSi[i, j, k]
    
    θ = DSi / (D + eps(0.0))

    D_linear_mortality, D_quadratic_mortality = mortality(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields)

    total_grazing = grazing(bgc.zooplankton, Val(:D), i, j, k, grid, bgc, clock, fields)

    return (total_grazing + D_linear_mortality + D_quadratic_mortality) * θ
end

@inline function calcite_production(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields)
    linear_mortality, quadratic_mortality = mortality(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields)

    total_grazing_loss = calcite_loss(bgc.zooplankton, Val(:P), i, j, k, grid, bgc, clock, fields)

    return total_grazing_loss + (linear_mortality + quadratic_mortality) / 2
end
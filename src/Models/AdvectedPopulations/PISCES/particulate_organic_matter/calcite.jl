@inline function (bgc::PISCES{<:Any, <:Any, <:Any, TwoCompartementCarbonIronParticles})(i, j, k, 
                                                                                        grid, 
                                                                                        val_name::Val{:CaCO₃},
                                                                                        clock, fields)
    R = rain_ratio(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields)
    
    phytoplankton_production = calcite_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields)

    return R * phytoplankton_production - dissolution
end

@inline function rain_ratio(pom::TwoCompartementCarbonIronParticles, i, j, k, grid, bgc, clock, fields)
    r = pom.base_rain_ratio

    T = @inbounds fields.T[i, j, k]
    PAR = @inbounds fields.PAR[i, j, k]
    zₘₓₗ = @inbounds fields.zₘₓₗ[i, j, k]

    L_CaCO₃ = coccolithophore_nutrient_limitation(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields)

    phytoplankton_concentration_factor = 
        coccolithophore_phytoplankton_factor(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields)

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

@inline function calcite_dissolution(poc::TwoCompartementCarbonIronParticles, i, j, k, grid, bgc, clock, fields)
    λ   = poc.base_calcite_dissolution_rate
    nca = poc.calcite_dissolution_exponent

    Ω     = @inbounds     fields.Ω[i, j, k]
    CaCO₃ = @inbounds fields.CaCO₃[i, j, k]

    ΔCaCO₃ = max(0, 1 - Ω)

    return λ * ΔCaCO₃ ^ nca * CaCO₃
end

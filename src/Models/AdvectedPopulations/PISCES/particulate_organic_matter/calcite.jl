@inline function (bgc::TwoCompartementPOCPISCES)(i, j, k, grid, val_name::Val{:CaCO₃}, clock, fields, auxiliary_fields)
    phytoplankton_production = calcite_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    dissolution = calcite_dissolution(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    
    return phytoplankton_production - dissolution
end

@inline function calcite_dissolution(poc::TwoCompartementCarbonIronParticles, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    λ   = poc.base_calcite_dissolution_rate
    nca = poc.calcite_dissolution_exponent

    Ω     = @inbounds auxiliary_fields.Ω[i, j, k]
    CaCO₃ = @inbounds       fields.CaCO₃[i, j, k]

    ΔCaCO₃ = max(0, 1 - Ω)

    return λ * ΔCaCO₃ ^ nca * CaCO₃
end

"""
    DissolvedOrganicCarbon

Parameterisation of dissolved organic matter which depends on a bacterial 
concentration.
"""
@kwdef struct DissolvedOrganicCarbon{FT, AP}
                             remineralisation_rate :: FT = 0.3/day   # 1 / s
             bacteria_concentration_depth_exponent :: FT = 0.684     # 
                  reference_bacteria_concentration :: FT = 1.0       # mmol C / m³
                           temperature_sensetivity :: FT = 1.066     #
# (1 / (mmol C / m³),  1 / (mmol C / m³),  1 / (mmol C / m³),  1 / (mmol C / m³) / s,  1 / (mmol C / m³) / s)
                            aggregation_parameters :: AP = (0.37, 102, 3530, 5095, 114) .* (10^-6 / day) 
end

required_biogeochemical_tracers(::DissolvedOrganicCarbon) = tuple(:DOC)

@inline function (bgc::PISCES{<:Any, <:Any, DissolvedOrganicCarbon})(i, j, k, grid, ::Val{:DOC}, clock, fields)
    phytoplankton_exudate = dissolved_exudate(bgc.phytoplankton, bgc, i, j, k, grid, bgc, clock, fields)
    upper_trophic_exudate = upper_trophic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)
    grazing_waste = organic_excretion(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)
    particulate_breakdown = degredation(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields)
    dissolved_breakdown = degredation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields)
    aggrgation_to_particles, = aggregation(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields)

    return (phytoplankton_exudate + upper_trophic_exudate + grazing_waste + particulate_breakdown 
            - dissolved_breakdown - aggregation_to_particles)
end

@inline function degredation(dom::DissolvedOrganicCarbon, i, j, k, grid, bgc, clock, fields)
    Bact_ref = dom.reference_bacteria_concentration
    b = dom.temperature_sensetivity
    λ = dom.remineralisation_rate

    T   = @inbounds   fields.T[i, j, k]
    DOC = @inbounds fields.DOC[i, j, k]

    f = b^T

    Bact = bacteria_concentration(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)

    LBact = bacteria_activity(bgc.zooplankton, i, j, k, grid, bgc, clock, fields)

    return λ * f * LBact * Bact / Bact_ref * DOC # differes from Aumont 2015 since the dimensions don't make sense 
end

@inline function aggregation(dom::DissolvedOrganicCarbon, i, j, k, grid, bgc, clock, fields)
    a₁, a₂, a₃, a₄, a₅ = dom.aggregation_parameters

    backgroound_shear = bgc.background_shear
    mixed_layer_shear = bgc.mixed_layer_shear

    z = znode(i, j, k, grid, Center(), Center(), Center())

    zₘₓₗ = @inbounds fields.zₘₓₗ[i, j, k]
    
    DOC = @inbounds fields.DOC[i, j, k]
    POC = @inbounds fields.POC[i, j, k]
    GOC = @inbounds fields.GOC[i, j, k]

    shear = ifelse(z < zₘₓₗ, backgroound_shear, mixed_layer_shear)
    
    Φ₁ = shear * (a₁ * DOC + a₂ * POC) * DOC
    Φ₂ = shear * (a₃ * GOC) * DOC
    Φ₃ = (a₄ * POC + a₅ * DOC) * DOC
    
    return Φ₁ + Φ₂ + Φ₃, Φ₁, Φ₂, Φ₃
end

@inline function aggregation_of_colloidal_iron(dom::DissolvedOrganicCarbon, i, j, k, grid, bgc, clock, fields)
    _, Φ₁, Φ₂, Φ₃ = aggregation(dom, i, j, k, grid, bgc, clock, fields)

    Fe′ = free_iron(bgc.iron, i, j, k, grid, bgc, clock, fields)
    ligand_iron = Fe - Fe′
    colloidal_iron = 0.5 * ligand_iron

    CgFe1 = (Φ₁ + Φ₃) * colloidal_iron / (DOC + eps(0.0))
    CgFe2 = Φ₂ * colloidal_iron / (DOC + eps(0.0))

    return CgFe1 + CgFe2, CgFe1, CgFe2
end

@inline function oxic_remineralisation(dom::DissolvedOrganicCarbon, i, j, k, grid, bgc, clock, fields)
    O₂ = @inbounds fields.O₂[i, j, k]

    ΔO₂ = anoxia_factor(bgc, O₂)

    degredation = degredation(dom, i, j, k, grid, bgc, clock, fields)

    return (1 - ΔO₂) * degredation
end

@inline function anoxic_remineralisation(dom::DissolvedOrganicCarbon, i, j, k, grid, bgc, clock, fields)
    O₂ = @inbounds fields.O₂[i, j, k]

    ΔO₂ = anoxia_factor(bgc, O₂)

    degredation = degredation(dom, i, j, k, grid, bgc, clock, fields)

    return ΔO₂ * degredation
end
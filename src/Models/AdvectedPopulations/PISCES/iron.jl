"""
    SimpleIron(; excess_scavenging_enhancement = 1000)

Parameterisation for iron evolution, not the "complex chemistry" model
of Aumount et al, 2015. Iron is scavenged (i.e. perminemtly removed from
the model) when the free iron concentration exeeds the ligand concentration
at a rate modified by `excess_scavenging_enhancement`.
"""
@kwdef struct SimpleIron{FT}
    excess_scavenging_enhancement :: FT = 1000 # unitless
end

@inline function (iron::SimpleIron)(::Val{:Fe}, bgc,
                                    x, y, z, t,
                                    P, D, Z, M, 
                                    PChl, DChl, PFe, DFe, DSi,
                                    DOC, POC, GOC, 
                                    SFe, BFe, PSi, 
                                    NO₃, NH₄, PO₄, Fe, Si, 
                                    CaCO₃, DIC, Alk, 
                                    O₂, T, S,
                                    zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    λ̄ = iron.excess_scavenging_enhancement

    λFe = iron_scavenging_rate(bgc.particulate_organic_matter, POC, GOC, CaCO₃, PSi)

    Fe′ = free_iron(iron, Fe, DOC, T)
    total_ligand_concentration = max(0.6, 0.09 * (DOC + 40) - 3)

    # terminal process which removes iron from the ocean
    ligand_aggregation = λ̄ * λFe * max(0, Fe - total_ligand_concentration) * Fe′

    # other aggregation
    colloidal_aggregation, = aggregation_of_colloidal_iron(iron, bgc.dissolved_organic_matter, bgc, z, DOC, POC, GOC, Fe, T, zₘₓₗ)

    aggregation = colloidal_aggregation + ligand_aggregation

    # scavening and bacterial uptake to particles
    scav = λFe * (POC + GOC) * Fe′

    BactFe = bacterial_iron_uptake(bgc.dissolved_organic_matter, z, Z, M, DOC, NO₃, NH₄, PO₄, Fe, T, zₘₓₗ, zₑᵤ)

    λPOC = specific_degredation_rate(bgc.particulate_organic_matter, bgc, O₂, T)

    # particle breakdown
    particulate_degredation = λPOC * SFe

    # consumption
    nanophytoplankton_consumption, = iron_uptake(bgc.nanophytoplankton, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T)
    diatom_consumption, = iron_uptake(bgc.diatoms, bgc, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T)

    consumption = nanophytoplankton_consumption + diatom_consumption

    # grazing waste - this is the excess non assimilated into zooplankton when they consume iron rich phytoplankton
    microzooplankton_waste = specific_non_assimilated_iron(bgc.microzooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, BFe, T) * Z
    mesozooplankton_waste  = specific_non_assimilated_iron(bgc.mesozooplankton, bgc, x, y, z, P, D, PFe, DFe, Z, POC, GOC, SFe, BFe, T) * M

    zooplankton_waste = microzooplankton_waste + mesozooplankton_waste

    # type in Aumount 2015, γ should not be present since DOC doesn't contain iron/there is no DOFe pool
    respiration_product = upper_trophic_respiration_product(bgc.mesozooplankton, M, T) * bgc.mesozooplankton.iron_ratio

    return zooplankton_waste + respiration_product + particulate_degredation - consumption - scav - aggregation - BactFe
end

@inline function free_iron(::SimpleIron, Fe, DOC, T)
    # maybe some of these numbers should be parameters
    ligands = max(0.6, 0.09 * (DOC + 40) - 3)
    K = exp(16.27 - 1565.7 / max(T + 273.15, 5))
    Δ = 1 + K * ligands - K * Fe

    return (-Δ + √(Δ^2 + 4K * Fe)) / 2K
end

# this should be dispatched on an abstract type if we implement complex chemistry
@inline function aggregation_of_colloidal_iron(iron::SimpleIron, dom, bgc, z, DOC, POC, GOC, Fe, T, zₘₓₗ)
    _, Φ₁, Φ₂, Φ₃ = aggregation(dom, bgc, z, DOC, POC, GOC, zₘₓₗ)

    Fe′ = free_iron(iron, Fe, DOC, T)
    ligand_iron = Fe - Fe′
    colloidal_iron = 0.5 * ligand_iron

    CgFe1 = (Φ₁ + Φ₃) * colloidal_iron / (DOC + eps(0.0))
    CgFe2 = Φ₂ * colloidal_iron / (DOC + eps(0.0))

    return CgFe1 + CgFe2, CgFe1, CgFe2
end

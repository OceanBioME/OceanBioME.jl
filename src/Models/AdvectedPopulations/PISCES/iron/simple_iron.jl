"""
    SimpleIron(; excess_scavenging_enhancement = 1000)

Parameterisation for iron evolution, not the "complex chemistry" model
of Aumount et al, 2015. Iron is scavenged (i.e. perminemtly removed from
the model) when the free iron concentration exeeds the ligand concentration
at a rate modified by `excess_scavenging_enhancement`.
"""
@kwdef struct SimpleIron{FT}
    excess_scavenging_enhancement :: FT = 1000.0 # unitless
     maximum_ligand_concentration :: FT = 0.6    # μmol Fe / m³
           dissolved_ligand_ratio :: FT = 0.09   # μmol Fe / mmol C
end

required_biogeochemical_tracers(::SimpleIron) = tuple(:Fe)

const SimpleIronPISCES = PISCES{<:Any, <:Any, <:Any, <:Any, <:Any, <:SimpleIron}

@inline function (bgc::SimpleIronPISCES)(i, j, k, grid, val_name::Val{:Fe}, clock, fields, auxiliary_fields)
    Fe = @inbounds fields.Fe[i, j, k]
    DOC = @inbounds fields.DOC[i, j, k]

    T = @inbounds fields.T[i, j, k]

    scavenging_rate = iron_scavenging_rate(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    free_iron_concentration = free_iron(bgc.iron, Fe, DOC, T)
    total_ligand_concentration = ligand_concentration(bgc.iron, DOC)

    ligand_aggregation_loss = ligand_aggregation(bgc.iron,
                                                 scavenging_rate,
                                                 Fe,
                                                 total_ligand_concentration,
                                                 free_iron_concentration)

    colloidal_aggregation, = aggregation_of_colloidal_iron(bgc.dissolved_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    scavenging = iron_scavenging(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    bacterial_uptake = bacterial_iron_uptake(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    small_particles = degredation(bgc.particulate_organic_matter, Val(:SFe), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    consumption = uptake(bgc.phytoplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    grazing_waste = non_assimilated_iron(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    upper_trophic_waste = upper_trophic_dissolved_iron(bgc.zooplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return iron_tendency(small_particles,
                         grazing_waste,
                         upper_trophic_waste,
                         consumption,
                         ligand_aggregation_loss,
                         colloidal_aggregation,
                         scavenging,
                         bacterial_uptake)
end

using Oceananigans.Grids: znode, Center

@kwdef struct MicroAndMeso{μ, M, FT}
                                             micro :: μ
                                              meso :: M

           microzooplankton_bacteria_concentration :: FT = 0.7
            mesozooplankton_bacteria_concentration :: FT = 1.4    
                    maximum_bacteria_concentration :: FT = 4.0       # mmol C / m³
             bacteria_concentration_depth_exponent :: FT = 0.684     # 

        doc_half_saturation_for_bacterial_activity :: FT = 417.0     # mmol C / m³
    nitrate_half_saturation_for_bacterial_activity :: FT = 0.03      # mmol N / m³
    ammonia_half_saturation_for_bacterial_activity :: FT = 0.003     # mmol N / m³
  phosphate_half_saturation_for_bacterial_activity :: FT = 0.003     # mmol P / m³
       iron_half_saturation_for_bacterial_activity :: FT = 0.01      # μmol Fe / m³
end

required_biogeochemical_tracers(zoo::MicroAndMeso) = (required_biogeochemical_tracers(zoo.micro, :Z)...,
                                                      required_biogeochemical_tracers(zoo.meso, :M)...)

@inline zooplankton_concentration(::Val{:Z}, i, j, k, fields) = @inbounds fields.Z[i, j, k]
@inline zooplankton_concentration(::Val{:M}, i, j, k, fields) = @inbounds fields.M[i, j, k]

@inline parameterisation(::Val{:Z}, zoo::MicroAndMeso) = zoo.micro
@inline parameterisation(::Val{:M}, zoo::MicroAndMeso) = zoo.meso

@inline predator_parameterisation(val_name, zoo) = nothing
@inline predator_parameterisation(::Val{:Z}, zoo::MicroAndMeso) = zoo.meso

@inline predator_name(val_name, zoo) = nothing
@inline predator_name(::Val{:Z}, zoo::MicroAndMeso) = Val(:M)

@inline grazing(::Nothing, ::Nothing, val_prey_name, i, j, k, grid, args...) = zero(grid)

@inline function (bgc::PISCES{<:Any, <:MicroAndMeso})(i, j, k, grid, val_name::Union{Val{:Z}, Val{:M}}, clock, fields, auxiliary_fields)
    zoo = parameterisation(val_name, bgc.zooplankton)

    net_production = growth_death(zoo, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    # M preying on Z
    predator_zoo = predator_parameterisation(val_name, bgc.zooplankton)
    val_predator_name = predator_name(val_name, bgc.zooplankton)

    predatory_grazing = grazing(predator_zoo, val_predator_name, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return net_production - predatory_grazing
end

@inline grazing(zoo::MicroAndMeso, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (grazing(zoo.micro, Val(:Z), val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + grazing(zoo.meso, Val(:M), val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline flux_feeding(zoo::MicroAndMeso, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
     (flux_feeding(zoo.micro, Val(:Z), val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
      + flux_feeding(zoo.meso, Val(:M), val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields))
 
@inline inorganic_excretion(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (inorganic_excretion(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + inorganic_excretion(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline organic_excretion(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (organic_excretion(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + organic_excretion(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline non_assimilated_iron(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (non_assimilated_iron(zoo.micro, Val(:Z), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + non_assimilated_iron(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline upper_trophic_excretion(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    upper_trophic_excretion(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline upper_trophic_respiration(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    upper_trophic_respiration(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     
@inline upper_trophic_dissolved_iron(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    upper_trophic_dissolved_iron(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    
@inline upper_trophic_fecal_production(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    upper_trophic_fecal_production(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline upper_trophic_fecal_iron_production(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    upper_trophic_fecal_iron_production(zoo.meso, Val(:M), i, j, k, grid, bgc, clock, fields, auxiliary_fields)

@inline function bacteria_concentration(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    bZ = zoo.microzooplankton_bacteria_concentration
    bM = zoo.mesozooplankton_bacteria_concentration
    a  = zoo.bacteria_concentration_depth_exponent

    z = znode(i, j, k, grid, Center(), Center(), Center())

    zₘₓₗ = @inbounds auxiliary_fields.zₘₓₗ[i, j, k]
    zₑᵤ  = @inbounds  auxiliary_fields.zₑᵤ[i, j, k]

    Z = @inbounds fields.Z[i, j, k]
    M = @inbounds fields.M[i, j, k]

    zₘ = min(zₘₓₗ, zₑᵤ)

    surface_bacteria = min(4, bZ * Z + bM * M)

    depth_factor = (zₘ / z) ^ a

    return ifelse(z >= zₘ, 1, depth_factor) * surface_bacteria
end

@inline function bacteria_activity(zoo::MicroAndMeso, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    K_DOC = zoo.doc_half_saturation_for_bacterial_activity
    K_NO₃ = zoo.nitrate_half_saturation_for_bacterial_activity
    K_NH₄ = zoo.ammonia_half_saturation_for_bacterial_activity
    K_PO₄ = zoo.phosphate_half_saturation_for_bacterial_activity
    K_Fe  = zoo.iron_half_saturation_for_bacterial_activity

    NH₄ = @inbounds fields.NH₄[i, j, k]
    NO₃ = @inbounds fields.NO₃[i, j, k]
    PO₄ = @inbounds fields.PO₄[i, j, k]
    Fe  = @inbounds  fields.Fe[i, j, k]
    DOC = @inbounds fields.DOC[i, j, k]

    DOC_limit = DOC / (DOC + K_DOC)

    L_N   = (K_NO₃ * NH₄ + K_NH₄ * NO₃) / (K_NO₃ * K_NH₄ + K_NO₃ * NH₄ + K_NH₄ * NO₃)

    L_PO₄ = PO₄ / (PO₄ + K_PO₄)

    L_Fe  = Fe / (Fe + K_Fe)

    # assuming typo in paper otherwise it doesn't make sense to formulate L_NH₄ like this
    limiting_quota = min(L_N, L_PO₄, L_Fe)

    return limiting_quota * DOC_limit
end

@inline calcite_loss(zoo::MicroAndMeso, val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (calcite_loss(zoo.micro, Val(:Z), val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + calcite_loss(zoo.meso, Val(:M), val_prey_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields))

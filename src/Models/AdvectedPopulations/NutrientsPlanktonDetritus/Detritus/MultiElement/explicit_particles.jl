# An explicit sinking version of MARBL's detritus model with dissolved/refractory/explicit particles racking carbon/nitrogen/phosphate

#####
##### Dissolved, refractory and particulate organic matter in carbon, nitrogen and phosphorus
#####
##### Semi-labile DOC/DON/DOP and refractory DOCr/DONr/DOPr, plus particulate POC/POP, particulate iron
##### PFe and biogenic silica bSi. Particulate nitrogen is diagnostic (PON = Q·POC — particulate N:C is
##### fixed, so it is not carried as a tracer). Detrital calcite is *not* here: it belongs to the
##### `ExplicitCalciumCarbonate` pool in the inorganic-carbon slot.
#####
##### `AbstractMultiElementRefractoryDissolvedDetritus` is the shared dissolved half — it does not care how particles
##### sink, so the six dissolved rate fields and their remineralisation functions are written once and
##### inherited by both the explicit-sinking type here and the implicit ballast one. The seam between them
##### is exactly the particulate remineralisation functions.
#####

abstract type AbstractMultiElementRefractoryDissolvedDetritus end

struct MultiElementRefractoryDissolvedParticulate{FT, R, F, SV} <: AbstractMultiElementRefractoryDissolvedDetritus
       semilabile_remineralisation_rate_light :: R   # 1/s, in lit water; NamedTuple over the elements
        semilabile_remineralisation_rate_dark :: R   # 1/s, in dark water
             refractory_remineralisation_rate :: R   # 1/s, refractory base rates

                        photodegradation_rate :: FT   # 1/s; surface UV augmentation
    particulate_organic_remineralisation_rate :: FT   # 1/s; POC/POP/PFe
     particulate_silica_remineralisation_rate :: FT   # 1/s; bSi

               refractory_production_fraction :: F    # share of dissolved production that is refractory
              particulate_refractory_fraction :: F    # share of particulate remineralisation → refractory

                           uv_reference_depth :: FT   # m; the depth the surface UV dose is spread over

                           sinking_velocities :: SV
end

const MultiElementRefractoryDissolvedDetritus_NPD =
    NutrientsPlanktonDetritus{<:Any, <:Nutrients, <:Any, <:AbstractMultiElementRefractoryDissolvedDetritus, <:Any}

const MultiElementRefractoryDissolvedParticulate_NPD =
    NutrientsPlanktonDetritus{<:Any, <:Nutrients, <:Any, <:MultiElementRefractoryDissolvedParticulate, <:Any}

required_biogeochemical_tracers(::MultiElementRefractoryDissolvedParticulate) =
    (:DOC, :DON, :DOP, :DOCr, :DONr, :DOPr, :POC, :POP, :bSi, :PFe)

required_biogeochemical_auxiliary_fields(::MultiElementRefractoryDissolvedParticulate) = tuple()  # PAR via the plankton

"""
    MultiElementRefractoryDissolvedParticulate(grid; kwargs...)

A multi-element detritus component for the `detritus` slot of a [`NutrientsPlanktonDetritus`](@ref)
model, carrying semi-labile dissolved organic carbon, nitrogen and phosphorus (`DOC`/`DON`/`DOP`),
their refractory counterparts (`DOCr`/`DONr`/`DOPr`), particulate organic carbon and phosphorus
(`POC`/`POP`), particulate iron (`PFe`) and biogenic silica (`bSi`), all sinking explicitly.

Particulate nitrogen is diagnostic (`PON = Q·POC`, particulate N:C being fixed) and detrital calcium
carbonate lives in the inorganic-carbon slot, so neither is a tracer here.

Semi-labile material remineralises at a light-dependent rate — a faster rate in lit water than in the
dark — and refractory material at a slow base rate augmented by photodegradation in the surface cell.
See [`MultiElementRefractoryDissolved`](@ref) for the implicit ballast-sinking sibling.

Keyword Arguments
=================

- `particulate_organic_remineralisation_rate_per_day`, `particulate_silica_remineralisation_rate_per_day`:
  first-order particulate remineralisation rates (1/day)
- `uv_reference_depth`: the depth over which the surface UV dose is spread (m)
- `sinking_speed`: the sinking speed of every particulate class (m/day)
- `open_bottom`: whether particulates sink out of the bottom of the domain
- see the constructor definition for the dissolved remineralisation rates and refractory fractions
"""
function MultiElementRefractoryDissolvedParticulate(grid::AbstractGrid{FT};
                                                    photodegradation_rate = 1/(18yr),
                                                    uv_reference_depth = 10.0,
                                                    light_rates = (carbon = 1/(15yr),    nitrogen = 1/(15yr),   phosphorus = 1/(60yr)),
                                                    dark_rates = (carbon = 1/(6yr),     nitrogen = 1/(5.5yr),  phosphorus = 1/(4.5yr)),
                                                    reference_rates = (carbon = 1/(16000yr), nitrogen = 1/(9500yr), phosphorus = 1/(5500yr)),
                                                    production_refractory_fraction = (carbon = 0.01,      nitrogen = 0.0115,      phosphorus = 0.003),
                                                    particulate_refractory_fraction = (carbon = 0.01*0.06, nitrogen = 0.0115*0.03, phosphorus = 0.003*0.06),
                                                    sinking_speed = 10.0,        # m/day
                                                    open_bottom = true) where FT

    sinking_velocities = setup_velocity_fields((; POC = convert(FT, sinking_speed / day)),
                                               grid, open_bottom; three_D = true).POC

    RT = NamedTuple{(:carbon, :nitrogen, :phosphorus), NTuple{3, FT}}

    SV = typeof(sinking_velocities)

    return MultiElementRefractoryDissolvedParticulate{FT, RT, RT, SV}(
        RT(light_rates),
        RT(dark_rates),
        RT(reference_rates),
        convert(FT, photodegradation_rate),
        convert(FT, particulate_organic_remineralisation_rate_per_day / day),
        convert(FT, particulate_silica_remineralisation_rate_per_day / day),
        RT(production_refractory_fraction),
        RT(particulate_refractory_fraction),
        convert(FT, uv_reference_depth),
        sinking_velocities
    )
end

#####
##### Remineralisation rates
#####
@subcolumn_average @inline light_dark_rate(PAR, rate_light, rate_dark) = 
    ifelse(PAR > one(PAR), rate_light, rate_dark)
    
@inline function refractory_rate(i, j, k, grid, aux, base_rate, photodegradation_rate, uv_reference_depth)
    I       = interface_par(i, j, k + 1, grid, aux) # if subcolumn info available returns subcolumn value
    Δz      = Δzᶜᶜᶜ(i, j, k, grid)
    surface = k == grid.Nz

    uv = uv_dependant_refractory_rate(I, surface, uv_reference_depth / Δz, photodegradation_rate)

    return base_rate + uv
end

@subcolumn_average @inline uv_dependant_refractory_rate(PAR, surface, uv_depth, rate) = 
    ifelse(surface & (PAR > one(PAR)),
           log(max(PAR, one(PAR))) * oftype(PAR, 0.4373) * uv_depth * rate,
           zero(PAR))

for (element, tracer, field) in ((:carbon,     :DOC, :doc),
                                 (:nitrogen,   :DON, :don),
                                 (:phosphorus, :DOP, :dop))
    @eval @inline function $(Symbol(field, :_remineralisation))(i, j, k, grid,
                                                                d::AbstractMultiElementRefractoryDissolvedDetritus,
                                                                bgc, fields, aux)
        @inbounds begin
            C = fields.$tracer[i, j, k]
            PAR = @preserve_subcolumns aux.PAR[i, j, k]
        end

        return C * light_dark_rate(PAR,
                                   d.semilabile_remineralisation_rate_light.$element,
                                   d.semilabile_remineralisation_rate_dark.$element)
    end

    @eval @inline function $(Symbol(field, :r_remineralisation))(i, j, k, grid,
                                                                 d::AbstractMultiElementRefractoryDissolvedDetritus,
                                                                 bgc, fields, aux)
        @inbounds C = fields.$(Symbol(tracer, :r))[i, j, k]

        return C * refractory_rate(i, j, k, grid, aux,
                                   d.refractory_remineralisation_rate.$element,
                                   d.photodegradation_rate, d.uv_reference_depth)
    end
end

@inline poc_remineralisation(i, j, k, grid, d::MultiElementRefractoryDissolvedParticulate, fields) =
    @inbounds fields.POC[i, j, k] * d.particulate_organic_remineralisation_rate

@inline pop_remineralisation(i, j, k, grid, d::MultiElementRefractoryDissolvedParticulate, fields) =
    @inbounds fields.POP[i, j, k] * d.particulate_organic_remineralisation_rate

@inline pfe_remineralisation(i, j, k, grid, d::MultiElementRefractoryDissolvedParticulate, fields) =
    @inbounds fields.PFe[i, j, k] * d.particulate_organic_remineralisation_rate

@inline bsi_remineralisation(i, j, k, grid, d::MultiElementRefractoryDissolvedParticulate, fields) =
    @inbounds fields.bSi[i, j, k] * d.particulate_silica_remineralisation_rate

@inline pon_remineralisation(i, j, k, grid, d::AbstractMultiElementRefractoryDissolvedDetritus, bgc, fields) =
    nitrogen_ratio(i, j, k, grid, bgc.plankton, bgc, fields) *
    poc_remineralisation(i, j, k, grid, d, fields)

@inline doc_production(i, j, k, grid, bgc, fields, aux) =
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux)

@inline don_production(i, j, k, grid, bgc, fields, aux) =
    nitrogen_ratio(i, j, k, grid, bgc.plankton, bgc, fields) *
    bgc.plankton.nitrogen_to_DON_fraction *
    doc_production(i, j, k, grid, bgc, fields, aux)

@inline poc_production(i, j, k, grid, bgc, fields, aux) =
    solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux)

@inline (bgc::MultiElementRefractoryDissolvedDetritus_NPD)(i, j, k, grid, ::Val{:DOC}, clock, fields, aux) =
    (1 - bgc.detritus.refractory_production_fraction.carbon) *
        doc_production(i, j, k, grid, bgc, fields, aux) -
    doc_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields, aux)

@inline (bgc::MultiElementRefractoryDissolvedDetritus_NPD)(i, j, k, grid, ::Val{:DOCr}, clock, fields, aux) =
    bgc.detritus.refractory_production_fraction.carbon *
        doc_production(i, j, k, grid, bgc, fields, aux) -
    docr_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields, aux) +
    bgc.detritus.particulate_refractory_fraction.carbon *
        poc_remineralisation(i, j, k, grid, bgc.detritus, fields)

@inline (bgc::MultiElementRefractoryDissolvedDetritus_NPD)(i, j, k, grid, ::Val{:DON}, clock, fields, aux) =
    (1 - bgc.detritus.refractory_production_fraction.nitrogen) *
        don_production(i, j, k, grid, bgc, fields, aux) -
    don_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields, aux)

@inline (bgc::MultiElementRefractoryDissolvedDetritus_NPD)(i, j, k, grid, ::Val{:DONr}, clock, fields, aux) =
    bgc.detritus.refractory_production_fraction.nitrogen *
        don_production(i, j, k, grid, bgc, fields, aux) -
    donr_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields, aux) +
    bgc.detritus.particulate_refractory_fraction.nitrogen *
        pon_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields)

@inline (bgc::MultiElementRefractoryDissolvedDetritus_NPD)(i, j, k, grid, ::Val{:DOP}, clock, fields, aux) =
    (1 - bgc.detritus.refractory_production_fraction.phosphorus) *
        dissolved_phosphorus_production(i, j, k, grid, bgc.plankton, bgc, fields, aux) -
    dop_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields, aux) -
    dissolved_phosphorus_uptake(i, j, k, grid, bgc.plankton, bgc, fields, aux) -
    dissolved_phosphorus_balance(i, j, k, grid, bgc.plankton, bgc, fields, aux)

@inline (bgc::MultiElementRefractoryDissolvedDetritus_NPD)(i, j, k, grid, ::Val{:DOPr}, clock, fields, aux) =
    bgc.detritus.refractory_production_fraction.phosphorus *
        dissolved_phosphorus_production(i, j, k, grid, bgc.plankton, bgc, fields, aux) -
    dopr_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields, aux) +
    bgc.detritus.particulate_refractory_fraction.phosphorus *
        pop_remineralisation(i, j, k, grid, bgc.detritus, fields)

@inline (bgc::MultiElementRefractoryDissolvedParticulate_NPD)(i, j, k, grid, ::Val{:POC}, clock, fields, aux) =
    poc_production(i, j, k, grid, bgc, fields, aux) -
    poc_remineralisation(i, j, k, grid, bgc.detritus, fields)

@inline (bgc::MultiElementRefractoryDissolvedParticulate_NPD)(i, j, k, grid, ::Val{:POP}, clock, fields, aux) =
    particulate_phosphorus_production(i, j, k, grid, bgc.plankton, bgc, fields, aux) -
    pop_remineralisation(i, j, k, grid, bgc.detritus, fields)

@inline (bgc::MultiElementRefractoryDissolvedParticulate_NPD)(i, j, k, grid, ::Val{:PFe}, clock, fields, aux) =
    particulate_iron_production(i, j, k, grid, bgc.plankton, bgc, fields, aux) -
    pfe_remineralisation(i, j, k, grid, bgc.detritus, fields)

@inline (bgc::MultiElementRefractoryDissolvedParticulate_NPD)(i, j, k, grid, ::Val{:bSi}, clock, fields, aux) =
    biogenic_silica_production(i, j, k, grid, bgc.plankton, bgc, fields, aux) -
    bsi_remineralisation(i, j, k, grid, bgc.detritus, fields)

# silicate is a nutrient, but biogenic silica lives here, so the detritus owns its tendency
@inline (bgc::MultiElementRefractoryDissolvedDetritus_NPD)(i, j, k, grid, ::Val{:Si}, clock, fields, aux) =
    silicon_dissolution(i, j, k, grid, bgc.plankton, bgc, fields, aux) +
    bsi_remineralisation(i, j, k, grid, bgc.detritus, fields)

#####
##### Remineralisation feeding back to the inorganic pools
#####
##### The refractory-diverted fractions go to DOCr/DONr/DOPr above, so what reaches the inorganic pools is
##### the remainder. Iron has no refractory pool.
#####

@inline inorganic_carbon_waste(i, j, k, grid, d::AbstractMultiElementRefractoryDissolvedDetritus, bgc, fields, aux) =
    doc_remineralisation(i, j, k, grid, d, bgc, fields, aux) +
    docr_remineralisation(i, j, k, grid, d, bgc, fields, aux) +
    (1 - d.particulate_refractory_fraction.carbon) *
        poc_remineralisation(i, j, k, grid, d, fields)

@inline inorganic_nitrogen_waste(i, j, k, grid, d::AbstractMultiElementRefractoryDissolvedDetritus, bgc, fields, aux) =
    don_remineralisation(i, j, k, grid, d, bgc, fields, aux) +
    donr_remineralisation(i, j, k, grid, d, bgc, fields, aux) +
    (1 - d.particulate_refractory_fraction.nitrogen) *
        pon_remineralisation(i, j, k, grid, d, bgc, fields)

@inline inorganic_phosphate_waste(i, j, k, grid, d::AbstractMultiElementRefractoryDissolvedDetritus, bgc, fields, aux) =
    dop_remineralisation(i, j, k, grid, d, bgc, fields, aux) +
    dopr_remineralisation(i, j, k, grid, d, bgc, fields, aux) +
    (1 - d.particulate_refractory_fraction.phosphorus) *
        pop_remineralisation(i, j, k, grid, d, fields)

@inline inorganic_iron_waste(i, j, k, grid, d::AbstractMultiElementRefractoryDissolvedDetritus, bgc, fields, aux) =
    pfe_remineralisation(i, j, k, grid, d, fields)

#####
##### Sinking
#####

for particle in (:POC, :POP, :PFe, :bSi)
    @eval biogeochemical_drift_velocity(bgc::MultiElementRefractoryDissolvedParticulate_NPD, ::Val{$(QuoteNode(particle))}) =
        bgc.detritus.sinking_velocities
end

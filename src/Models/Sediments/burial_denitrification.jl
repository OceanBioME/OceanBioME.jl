#####
##### Burial, sedimentary denitrification and anaerobic remineralisation
#####
##### A sediment which closes the bottom boundary for a multi-element detritus. It is fed the sinking flux
##### of particulate organic carbon and phosphorus, biogenic silica, calcite and particulate iron reaching
##### the sea floor, buries part of each, and returns the rest to the bottom cell through the routes the
##### water column uses: organic carbon → dissolved inorganic carbon and refractory dissolved organic
##### carbon, organic nitrogen → ammonium and refractory dissolved organic nitrogen, organic phosphorus →
##### phosphate and its refractory pool, biogenic silica → silicate, calcite → inorganic carbon and
##### alkalinity. Particulate iron is removed entirely.
#####
##### Every method returns an AREA flux (mmol m⁻² s⁻¹); the framework divides by the bottom-cell thickness,
##### so the expressions here are independent of the vertical grid.
#####

# `IS` (implicit sinking) selects where the floor fluxes come from. The sediment is constructed
# independently of the biogeochemistry, so it cannot see the detritus type and has to be told:
#   IS = false — the particulates are tracers and the framework computes their floor flux for us
#   IS = true  — there are no particulate tracers; the ballast sweep already has the exact floor flux and
#                publishes it as an auxiliary field, which we read directly
struct BurialDenitrification{FT, F, IS} <: AbstractSedimentBiogeochemistry
    # particulate organic burial efficiency, Dunne et al. (2007)
             particulate_burial_intercept :: FT
           particulate_burial_coefficient :: FT
       particulate_burial_half_saturation :: FT  # mmol m⁻² d⁻¹
        organic_carbon_burial_coefficient :: FT
      organic_nitrogen_burial_coefficient :: FT
    organic_phosphorus_burial_coefficient :: FT
          maximum_organic_burial_fraction :: FT

    # biogenic silica burial efficiency, Ragueneau et al. (2000)
        silica_burial_flux_threshold :: FT  # mmol m⁻² d⁻¹
    silica_burial_fraction_high_flux :: FT
     silica_burial_fraction_low_flux :: FT
           silica_burial_coefficient :: FT
      maximum_silica_burial_fraction :: FT

    # calcite is preserved where the bottom water is above this saturation state
    calcite_burial_saturation_threshold :: FT

    # sedimentary denitrification, Bohlen et al. (2012)
    sedimentary_denitrification_coefficient :: FT
      sedimentary_denitrification_intercept :: FT
          sedimentary_denitrification_scale :: FT
           sedimentary_denitrification_base :: FT

    # the combined denitrification rate is capped at the nitrate present over this horizon
    nitrate_depletion_timescale :: FT  # s

    # anaerobic remineralisation by oxidants other than oxygen and nitrate, Soetaert et al. (1996)
    other_remineralisation_intercept :: FT
      other_remineralisation_maximum :: FT
          anoxic_bottom_water_oxygen :: FT  # mmol/m³, below which ALL remineralisation is anaerobic

    # shared with the water column — these MUST match the plankton, detritus and oxygen components
                    nitrogen_to_carbon :: FT
       particulate_refractory_fraction :: F
    denitrification_carbon_to_nitrogen :: FT
     remineralisation_carbon_to_oxygen :: FT
                        oxygen_minimum :: FT
                  oxygen_minimum_range :: FT

    # the sea-floor organic remineralisation drives the ligand exactly as the water column does, so these
    # must match the iron cycle
    ligand_production_ratio :: FT
    ligand_degradation_rate :: FT
end

"""
    BurialDenitrification([FT = Float64;] kwargs...)
    BurialDenitrificationSediment(grid; kwargs...)

A sediment biogeochemistry which buries part of the organic and mineral flux reaching the sea floor and
returns the rest to the bottom cell, with sedimentary denitrification and anaerobic remineralisation.
Designed to close the bottom boundary of a model using a multi-element detritus such as
[`MultiElementRefractoryDissolvedParticulate`](@ref).

Organic burial efficiency follows Dunne et al. (2007), biogenic silica follows Ragueneau et al. (2000),
sedimentary denitrification follows Bohlen et al. (2012), and remineralisation by oxidants other than
oxygen and nitrate follows Soetaert et al. (1996). Calcite is buried only where the overlying water is
above the calcite saturation threshold, and particulate iron is removed entirely.

Use `BurialDenitrificationSediment(grid; kwargs...)` to build the coupled sediment model.

Keyword Arguments
=================

- `particulate_burial_*`, `organic_*_burial_coefficient`, `maximum_organic_burial_fraction`: the organic
  burial efficiency
- `silica_burial_*`, `maximum_silica_burial_fraction`: the biogenic silica burial efficiency
- `calcite_burial_saturation_threshold`: the saturation state above which calcite is preserved
- `sedimentary_denitrification_*`, `nitrate_depletion_timescale`: sedimentary denitrification and the
  cap that stops it exhausting the bottom-cell nitrate
- `other_remineralisation_*`, `anoxic_bottom_water_oxygen`: anaerobic remineralisation
- `nitrogen_to_carbon`, `particulate_refractory_fraction`, `denitrification_carbon_to_nitrogen`,
  `remineralisation_carbon_to_oxygen`, `oxygen_minimum`, `oxygen_minimum_range`,
  `ligand_production_ratio`, `ligand_degradation_rate`: shared with the water column, and **must** match
  the values used by the other components
- `implicit_sinking`: `true` when the detritus solves its particulates as a diagnostic flux rather than
  carrying them as tracers
"""
function BurialDenitrification(FT = Float64;
                               particulate_burial_intercept          = 0.013,
                               particulate_burial_coefficient        = 0.53,
                               particulate_burial_half_saturation    = 7.0,    # mmol m⁻² d⁻¹
                               organic_carbon_burial_coefficient     = 2.54,
                               organic_nitrogen_burial_coefficient   = 0.5,
                               organic_phosphorus_burial_coefficient = 0.36,
                               maximum_organic_burial_fraction       = 0.8,

                               silica_burial_flux_threshold     = 2.0,   # mmol m⁻² d⁻¹
                               silica_burial_fraction_high_flux = 0.2,
                               silica_burial_fraction_low_flux  = 0.04,
                               silica_burial_coefficient        = 1.53,
                               maximum_silica_burial_fraction   = 1.0,

                               calcite_burial_saturation_threshold = 0.89,

                               sedimentary_denitrification_coefficient = 1.0,
                               sedimentary_denitrification_intercept   = 0.06,
                               sedimentary_denitrification_scale       = 0.19,
                               sedimentary_denitrification_base        = 0.99,

                               nitrate_depletion_timescale = 10day,

                               other_remineralisation_intercept = 0.1,
                               other_remineralisation_maximum   = 0.5,
                               anoxic_bottom_water_oxygen       = 1.0,   # mmol/m³

                               nitrogen_to_carbon = 16/117,
                               particulate_refractory_fraction =
                                   (carbon = 0.01 * 0.06, nitrogen = 0.0115 * 0.03, phosphorus = 0.003 * 0.06),
                               denitrification_carbon_to_nitrogen = 117/136,
                               remineralisation_carbon_to_oxygen  = 117/138,
                               oxygen_minimum                     = 5.0,   # mmol/m³
                               oxygen_minimum_range               = 5.0,   # mmol/m³

                               ligand_production_ratio = 1.0e-4,
                               ligand_degradation_rate = 9.4e-5,

                               implicit_sinking = false)

    refractory = map(x -> convert(FT, x), particulate_refractory_fraction)

    return BurialDenitrification{FT, typeof(refractory), implicit_sinking}(
        convert(FT, particulate_burial_intercept),
        convert(FT, particulate_burial_coefficient),
        convert(FT, particulate_burial_half_saturation),
        convert(FT, organic_carbon_burial_coefficient),
        convert(FT, organic_nitrogen_burial_coefficient),
        convert(FT, organic_phosphorus_burial_coefficient),
        convert(FT, maximum_organic_burial_fraction),
        convert(FT, silica_burial_flux_threshold),
        convert(FT, silica_burial_fraction_high_flux),
        convert(FT, silica_burial_fraction_low_flux),
        convert(FT, silica_burial_coefficient),
        convert(FT, maximum_silica_burial_fraction),
        convert(FT, calcite_burial_saturation_threshold),
        convert(FT, sedimentary_denitrification_coefficient),
        convert(FT, sedimentary_denitrification_intercept),
        convert(FT, sedimentary_denitrification_scale),
        convert(FT, sedimentary_denitrification_base),
        convert(FT, nitrate_depletion_timescale),
        convert(FT, other_remineralisation_intercept),
        convert(FT, other_remineralisation_maximum),
        convert(FT, anoxic_bottom_water_oxygen),
        convert(FT, nitrogen_to_carbon),
        refractory,
        convert(FT, denitrification_carbon_to_nitrogen),
        convert(FT, remineralisation_carbon_to_oxygen),
        convert(FT, oxygen_minimum),
        convert(FT, oxygen_minimum_range),
        convert(FT, ligand_production_ratio),
        convert(FT, ligand_degradation_rate))
end

BurialDenitrificationSediment(grid; FT = eltype(grid), kwargs...) =
    BiogeochemicalSediment(grid, BurialDenitrification(FT; kwargs...))

# the prognostic fields are burial and nitrogen-loss inventories (mmol m⁻²), so their tendencies are the
# burial fluxes and the column element budget closes. The buried mass leaves through the open bottom.
required_sediment_fields(::BurialDenitrification) =
    (:buried_C, :buried_N, :buried_P, :buried_Si, :buried_Fe, :buried_CaCO₃, :denitrified_N)

coupled_tracers(::BurialDenitrification) =
    (:NO₃, :NH₄, :PO₄, :Si, :DIC, :Alk, :O₂, :DOCr, :DONr, :DOPr, :Lig)

# explicit sinking — the framework computes the floor flux from the particulate tracers
required_tracers(::BurialDenitrification{<:Any, <:Any, false}) = (:O₂, :NO₃, :Ω, :O₂_consumption_scale)
sinking_fluxes(::BurialDenitrification{<:Any, <:Any, false})   = (:POC, :POP, :bSi, :CaCO₃, :PFe)

# implicit sinking — there are no particulate tracers; the sweep's floor flux is already published
required_tracers(::BurialDenitrification{<:Any, <:Any, true}) =
    (:O₂, :NO₃, :Ω, :O₂_consumption_scale, :POC_floor, :POP_floor, :bSi_floor, :CaCO₃_floor, :PFe_floor)
sinking_fluxes(::BurialDenitrification{<:Any, <:Any, true}) = ()

# ...so the tendencies below are written once against `tracked.POC` and the floor fluxes are remapped to
# those names once, here, at the call boundary
@inline particulate_floor_fluxes(::BurialDenitrification{<:Any, <:Any, false}, tracked) = tracked

@inline particulate_floor_fluxes(::BurialDenitrification{<:Any, <:Any, true}, t) =
    (POC = t.POC_floor, POP = t.POP_floor, bSi = t.bSi_floor, CaCO₃ = t.CaCO₃_floor, PFe = t.PFe_floor,
     O₂ = t.O₂, NO₃ = t.NO₃, Ω = t.Ω, O₂_consumption_scale = t.O₂_consumption_scale)

@inline (s::BurialDenitrification)(i, j, grid, val_name::Val, clock, fields, tracked) =
    _sediment_tendency(s, i, j, grid, val_name, clock, fields, particulate_floor_fluxes(s, tracked))

#####
##### Burial fractions
#####
##### Everything done to a pool is `F·f` buried and `F·(1 − f)` returned, so the fraction is the only thing
##### worth naming. Each takes an already-clamped floor flux: every rate vanishes at zero flux, so clamping
##### once in the tendency reproduces the guard branchlessly and keeps the returned remineralisation ≥ 0.
#####

@inline floor_flux(F) = max(F, zero(F))

@inline burial_flux_per_day(F) = F * day                     # mmol m⁻² d⁻¹
@inline remin_flux_per_year(F) = F * 1e-4 * 365day           # mmol cm⁻² yr⁻¹

@inline particulate_burial_fraction(a, b, k, F̂) = a + b * F̂^2 / (k + F̂)^2

@inline silica_burial_fraction(f_high, f_low, F̂_thres, F̂) = ifelse(F̂ > F̂_thres, f_high, f_low)

@inline particulate_burial_fraction(s::BurialDenitrification, F_POC) =
    particulate_burial_fraction(s.particulate_burial_intercept, s.particulate_burial_coefficient,
                                s.particulate_burial_half_saturation, burial_flux_per_day(F_POC))

@inline organic_carbon_burial_fraction(s::BurialDenitrification, F_POC) =
    min(s.maximum_organic_burial_fraction,
        s.organic_carbon_burial_coefficient * particulate_burial_fraction(s, F_POC))

# nitrogen burial is not its own efficiency: it is a coefficient times the buried carbon's nitrogen. The
# particulate N:C is fixed, so as a fraction of the nitrogen floor flux it is just that coefficient times
# the carbon fraction.
@inline organic_nitrogen_burial_fraction(s::BurialDenitrification, F_POC) =
    s.organic_nitrogen_burial_coefficient * organic_carbon_burial_fraction(s, F_POC)

# phosphorus reuses the carbon burial efficiency with its own coefficient and the same cap
@inline organic_phosphorus_burial_fraction(s::BurialDenitrification, F_POC) =
    min(s.maximum_organic_burial_fraction,
        s.organic_phosphorus_burial_coefficient * particulate_burial_fraction(s, F_POC))

@inline biogenic_silica_burial_fraction(s::BurialDenitrification, F_bSi) =
    min(s.maximum_silica_burial_fraction,
        s.silica_burial_coefficient * silica_burial_fraction(s.silica_burial_fraction_high_flux,
                                                             s.silica_burial_fraction_low_flux,
                                                             s.silica_burial_flux_threshold,
                                                             burial_flux_per_day(F_bSi)))

@inline calcite_burial_fraction(s::BurialDenitrification, Ω) = Ω > s.calcite_burial_saturation_threshold

# all the particulate iron reaching the floor is removed
@inline iron_burial_fraction(::BurialDenitrification, F_PFe) = one(F_PFe)

#####
##### RedoxOxygen pathways (area fluxes)
#####

@inline function sedimentary_denitrification(s::BurialDenitrification, F_POC, O₂, NO₃)
    a = s.sedimentary_denitrification_intercept
    b = s.sedimentary_denitrification_scale
    r = s.sedimentary_denitrification_base

    return s.sedimentary_denitrification_coefficient * F_POC * (a + b * r^(O₂ - NO₃))
end

# below `anoxic_bottom_water_oxygen` ALL non-buried remineralisation is denitrification plus other oxidants
@inline function other_remineralisation(s::BurialDenitrification, F_POC, O₂, NO₃)
    F̃ = remin_flux_per_year(F_POC)

    available   = F_POC * (1 - organic_carbon_burial_fraction(s, F_POC))
    denitrified = sedimentary_denitrification(s, F_POC, O₂, NO₃) * s.denitrification_carbon_to_nitrogen
    anoxic      = available - denitrified

    soetaert = min(min(s.other_remineralisation_intercept + F̃, s.other_remineralisation_maximum) * available,
                   anoxic)

    return ifelse(O₂ < s.anoxic_bottom_water_oxygen, anoxic, soetaert)
end

# the bottom cell's water-column denitrification responds to the sediment: the returned organic carbon
# raises it, and the sediment's own anaerobic pathways are subtracted. The suboxic fraction reads the same
# bottom-cell oxygen the interior kernel used, so this delta is exact.
@inline function denitrification_response(s::BurialDenitrification, F_POC, O₂, NO₃)
    ω = suboxic_fraction(O₂, s.oxygen_minimum, s.oxygen_minimum_range)

    oxidised = (1 - s.particulate_refractory_fraction.carbon) *
               F_POC * (1 - organic_carbon_burial_fraction(s, F_POC))

    return ω * ((oxidised - other_remineralisation(s, F_POC, O₂, NO₃)) / s.denitrification_carbon_to_nitrogen -
                sedimentary_denitrification(s, F_POC, O₂, NO₃))
end

@inline nitrogen_loss(s::BurialDenitrification, F_POC, O₂, NO₃) =
    sedimentary_denitrification(s, F_POC, O₂, NO₃) + denitrification_response(s, F_POC, O₂, NO₃)

# the cap that stops denitrification exhausting the bottom-cell nitrate over the depletion horizon. The
# rate is an area flux, so it is divided by the bottom-cell thickness — the coupling's own divisor — to
# recover the concentration consumed.
@inline function denitrification_clamp(s::BurialDenitrification, F_POC, O₂, NO₃, Δz)
    rate = s.nitrate_depletion_timescale * nitrogen_loss(s, F_POC, O₂, NO₃) / Δz

    NO₃⁺ = max(NO₃, zero(NO₃))

    return ifelse(rate > NO₃⁺, NO₃⁺ / rate, one(NO₃⁺))   # divides only when rate > NO₃⁺ ≥ 0, so never 0/0
end

# the bottom-cell thickness the cap needs. On an immersed grid the floor is not always the first cell, so
# walk up to the first live one — the same search the framework uses to gather the tracked fields, so this
# is the identical cell. A plain grid skips nothing. Tests calling a tendency directly pass `nothing`,
# meaning a single box with no vertical structure, where the cap is inert.
@inline function bottom_Δz(grid::ImmersedBoundaryGrid, i, j)
    Nz = size(grid, 3)
    k = 1

    while immersed_cell(i, j, k, grid.underlying_grid, grid.immersed_boundary) && k < Nz
        k += 1
    end

    return Δzᶜᶜᶠ(i, j, k, grid)
end

@inline bottom_Δz(grid, i, j) = Δzᶜᶜᶠ(i, j, 1, grid)
@inline bottom_Δz(::Nothing, i, j) = Inf

# the returned organic carbon is oxidised aerobically except the part respired on nitrate and by other
# oxidants. The capped denitrification is used, since the cap is applied before the oxygen budget, while
# the other-oxidant term was computed before it.
@inline function oxygen_consumption(s::BurialDenitrification, F_POC, O₂, NO₃, work)
    η = aerobic_fraction(O₂, s.oxygen_minimum, s.oxygen_minimum_range)

    oxidised = (1 - s.particulate_refractory_fraction.carbon) *
               F_POC * (1 - organic_carbon_burial_fraction(s, F_POC)) -
               work * sedimentary_denitrification(s, F_POC, O₂, NO₃) * s.denitrification_carbon_to_nitrogen -
               other_remineralisation(s, F_POC, O₂, NO₃)

    return η * oxidised / s.remineralisation_carbon_to_oxygen
end

@inline oxygen_consumption(s::BurialDenitrification, F_POC, O₂, NO₃) =
    oxygen_consumption(s, F_POC, O₂, NO₃, one(F_POC))

#####
##### Tendencies
#####

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:buried_C}, clock, fields, tracked)
    @inbounds F_POC = floor_flux(tracked.POC[i, j, 1])

    return F_POC * organic_carbon_burial_fraction(s, F_POC)
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:buried_N}, clock, fields, tracked)
    @inbounds F_POC = floor_flux(tracked.POC[i, j, 1])

    return s.nitrogen_to_carbon * F_POC * organic_nitrogen_burial_fraction(s, F_POC)
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:buried_P}, clock, fields, tracked)
    @inbounds begin
        F_POC = floor_flux(tracked.POC[i, j, 1])
        F_POP = floor_flux(tracked.POP[i, j, 1])
    end

    return F_POP * organic_phosphorus_burial_fraction(s, F_POC)
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:buried_Si}, clock, fields, tracked)
    @inbounds F_bSi = floor_flux(tracked.bSi[i, j, 1])

    return F_bSi * biogenic_silica_burial_fraction(s, F_bSi)
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:buried_Fe}, clock, fields, tracked)
    @inbounds F_PFe = floor_flux(tracked.PFe[i, j, 1])

    return F_PFe * iron_burial_fraction(s, F_PFe)
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:buried_CaCO₃}, clock, fields, tracked)
    @inbounds begin
        F_CaCO₃ = floor_flux(tracked.CaCO₃[i, j, 1])
        Ω       = tracked.Ω[i, j, 1]
    end

    return F_CaCO₃ * calcite_burial_fraction(s, Ω)
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:denitrified_N}, clock, fields, tracked)
    @inbounds begin
        F_POC = floor_flux(tracked.POC[i, j, 1])
        O₂    = tracked.O₂[i, j, 1]
        NO₃   = tracked.NO₃[i, j, 1]
    end

    work = denitrification_clamp(s, F_POC, O₂, NO₃, bottom_Δz(grid, i, j))

    return work * nitrogen_loss(s, F_POC, O₂, NO₃)
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:DIC}, clock, fields, tracked)
    @inbounds begin
        F_POC   = floor_flux(tracked.POC[i, j, 1])
        F_CaCO₃ = floor_flux(tracked.CaCO₃[i, j, 1])
        Ω       = tracked.Ω[i, j, 1]
    end

    return (1 - s.particulate_refractory_fraction.carbon) *
               F_POC * (1 - organic_carbon_burial_fraction(s, F_POC)) +
           F_CaCO₃ * (1 - calcite_burial_fraction(s, Ω))
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:DOCr}, clock, fields, tracked)
    @inbounds F_POC = floor_flux(tracked.POC[i, j, 1])

    return s.particulate_refractory_fraction.carbon * F_POC * (1 - organic_carbon_burial_fraction(s, F_POC))
end

# the sea-floor organic remineralisation drives the ligand exactly as the water-column one does, so the net
# benthic ligand is the production ratio less the degradation rate, times the returned organic carbon
@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:Lig}, clock, fields, tracked)
    @inbounds F_POC = floor_flux(tracked.POC[i, j, 1])

    POC_remin = F_POC * (1 - organic_carbon_burial_fraction(s, F_POC))

    return (s.ligand_production_ratio - s.ligand_degradation_rate) * POC_remin
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:NH₄}, clock, fields, tracked)
    @inbounds F_POC = floor_flux(tracked.POC[i, j, 1])

    F_PON = s.nitrogen_to_carbon * F_POC

    return (1 - s.particulate_refractory_fraction.nitrogen) *
           F_PON * (1 - organic_nitrogen_burial_fraction(s, F_POC))
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:DONr}, clock, fields, tracked)
    @inbounds F_POC = floor_flux(tracked.POC[i, j, 1])

    F_PON = s.nitrogen_to_carbon * F_POC

    return s.particulate_refractory_fraction.nitrogen *
           F_PON * (1 - organic_nitrogen_burial_fraction(s, F_POC))
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:PO₄}, clock, fields, tracked)
    @inbounds begin
        F_POC = floor_flux(tracked.POC[i, j, 1])
        F_POP = floor_flux(tracked.POP[i, j, 1])
    end

    return (1 - s.particulate_refractory_fraction.phosphorus) *
           F_POP * (1 - organic_phosphorus_burial_fraction(s, F_POC))
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:DOPr}, clock, fields, tracked)
    @inbounds begin
        F_POC = floor_flux(tracked.POC[i, j, 1])
        F_POP = floor_flux(tracked.POP[i, j, 1])
    end

    return s.particulate_refractory_fraction.phosphorus *
           F_POP * (1 - organic_phosphorus_burial_fraction(s, F_POC))
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:Si}, clock, fields, tracked)
    @inbounds F_bSi = floor_flux(tracked.bSi[i, j, 1])

    return F_bSi * (1 - biogenic_silica_burial_fraction(s, F_bSi))
end

@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:NO₃}, clock, fields, tracked)
    @inbounds begin
        F_POC = floor_flux(tracked.POC[i, j, 1])
        O₂    = tracked.O₂[i, j, 1]
        NO₃   = tracked.NO₃[i, j, 1]
    end

    work = denitrification_clamp(s, F_POC, O₂, NO₃, bottom_Δz(grid, i, j))

    return -work * nitrogen_loss(s, F_POC, O₂, NO₃)
end

# alkalinity is assembled inside the pointwise kernel, before the sediment fluxes are added, so the
# sediment must supply its own alkalinity flux
@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:Alk}, clock, fields, tracked)
    @inbounds begin
        F_POC   = floor_flux(tracked.POC[i, j, 1])
        F_CaCO₃ = floor_flux(tracked.CaCO₃[i, j, 1])
        O₂      = tracked.O₂[i, j, 1]
        NO₃     = tracked.NO₃[i, j, 1]
        Ω       = tracked.Ω[i, j, 1]
    end

    F_PON = s.nitrogen_to_carbon * F_POC
    work  = denitrification_clamp(s, F_POC, O₂, NO₃, bottom_Δz(grid, i, j))

    return work * nitrogen_loss(s, F_POC, O₂, NO₃) +
           (1 - s.particulate_refractory_fraction.nitrogen) *
               F_PON * (1 - organic_nitrogen_burial_fraction(s, F_POC)) +
           2 * F_CaCO₃ * (1 - calcite_burial_fraction(s, Ω))
end

# the benthic oxygen consumption is scaled by the same factor as the water-column consumption; there is no
# production at the floor to scale
@inline function _sediment_tendency(s::BurialDenitrification, i, j, grid, ::Val{:O₂}, clock, fields, tracked)
    @inbounds begin
        F_POC = floor_flux(tracked.POC[i, j, 1])
        O₂    = tracked.O₂[i, j, 1]
        NO₃   = tracked.NO₃[i, j, 1]
        scale = tracked.O₂_consumption_scale[i, j, 1]
    end

    work = denitrification_clamp(s, F_POC, O₂, NO₃, bottom_Δz(grid, i, j))

    return -scale * oxygen_consumption(s, F_POC, O₂, NO₃, work)
end

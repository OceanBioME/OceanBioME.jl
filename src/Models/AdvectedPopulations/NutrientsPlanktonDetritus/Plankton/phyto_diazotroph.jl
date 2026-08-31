using Oceananigans.Fields: ZeroField

#####
##### A single photoautotrophic functional group
#####

"""
    Photoautotroph([FT = Float64;] kwargs...)

The growth and loss parameters of a single photoautotrophic functional group. Groups are the
building blocks of the [`PhytoDiazotroph`](@ref) plankton component and are not a `plankton`
component themselves.

Growth is the product of a maximum rate, a light limitation, and the Liebig (minimum) nutrient
limitation over the group's `nutrient_half_saturations` — whose keys therefore set which nutrients
limit this group. A group with no `nitrate` (or `ammonia`) key is not nitrogen limited, which is how
diazotrophs are distinguished from other phytoplankton.

Keyword Arguments
=================

- `nutrient_half_saturations`: a `NamedTuple` of half-saturation constants keyed by limiting nutrient
  (`nitrate`, `ammonia`, `phosphate`, and/or `iron`); its keys set which nutrients limit growth
  (mmol X / m³)
- `light_half_saturation`: the light half-saturation (W / m²)
- `light_limitation`: the light-limitation functional form, e.g. [`ExponentialLightLimitation`](@ref)
- `maximum_growth_rate`: the light- and nutrient-replete gross growth rate (1 / s)
- `mortality_rate`: the linear mortality rate (1 / s)
- `quadratic_mortality_rate`: the density dependent mortality rate, zero by default, which closes the
  model where grazers are not resolved (1 / s / mmol N / m³)
- `iron_ratio`: the iron content of the group's biomass (mol Fe / mol N)
- `chlorophyll_ratio`: the chlorophyll content of the group's biomass (g Chl / mol N)
"""
struct Photoautotroph{LN, TT, FT, LL} <: AbstractPlankton{LN}
    nutrient_half_saturations :: NamedTuple{LN, TT}

        light_half_saturation :: FT
             light_limitation :: LL

          maximum_growth_rate :: FT
               mortality_rate :: FT
     quadratic_mortality_rate :: FT

                   iron_ratio :: FT
            chlorophyll_ratio :: FT
end

function Photoautotroph(FT = Float64;
                        nutrient_half_saturations = (nitrate = 0.1,      # mmol N / m³
                                                     phosphate = 0.03,   # mmol P / m³
                                                     iron = 8e-5),       # mmol Fe / m³
                        light_half_saturation = 2.17,                    # W / m²
                        light_limitation = ExponentialLightLimitation(),
                        maximum_growth_rate = 0.8 / day,                 # 1 / s
                        mortality_rate = 0.06 / day,                     # 1 / s
                        quadratic_mortality_rate = 0,                    # 1 / s / mmol N / m³
                        iron_ratio = 3.3125e-5,                          # mol Fe / mol N
                        chlorophyll_ratio = 1.137)                       # g Chl / mol N

    limiting_nutrients = keys(nutrient_half_saturations)

    nutrient_half_saturations = NamedTuple{limiting_nutrients}(map(val->convert(FT, val), values(nutrient_half_saturations)))

    return Photoautotroph(nutrient_half_saturations,
                          convert(FT, light_half_saturation),
                          light_limitation,
                          convert(FT, maximum_growth_rate),
                          convert(FT, mortality_rate),
                          convert(FT, quadratic_mortality_rate),
                          convert(FT, iron_ratio),
                          convert(FT, chlorophyll_ratio))
end

@inline limiting_nutrients(::Photoautotroph{LN}) where LN = LN

@inline nutrient_half_saturations(group::Photoautotroph, ::Val{:N})   = group.nutrient_half_saturations.nitrate
@inline nutrient_half_saturations(group::Photoautotroph, ::Val{:NO₃}) = group.nutrient_half_saturations.nitrate
@inline nutrient_half_saturations(group::Photoautotroph, ::Val{:NH₄}) = group.nutrient_half_saturations.ammonia
@inline nutrient_half_saturations(group::Photoautotroph, ::Val{:PO₄}) = group.nutrient_half_saturations.phosphate
@inline nutrient_half_saturations(group::Photoautotroph, ::Val{:Fe})  = group.nutrient_half_saturations.iron

"""Light and nutrient limited gross specific growth rate (1 / s) of `group`."""
@inline function specific_growth_rate(i, j, k, grid, group::Photoautotroph, bgc, fields, auxiliary_fields)
    μ₀ = group.maximum_growth_rate
    kₗ = group.light_half_saturation

    PAR = @inbounds auxiliary_fields.PAR[i, j, k]

    Lₙ = nutrient_limitation(i, j, k, grid, bgc.nutrients, group, bgc, fields, auxiliary_fields)
    Lₗ = light_limitation(group.light_limitation, PAR, kₗ)

    return μ₀ * Lₗ * Lₙ
end

"""Linear plus density dependent mortality (mmol N / m³ / s) of a `group` with biomass `B`."""
@inline mortality(group::Photoautotroph, B) = (group.mortality_rate + group.quadratic_mortality_rate * B) * B

Adapt.adapt_structure(to, group::Photoautotroph) =
    Photoautotroph(adapt(to, group.nutrient_half_saturations),
                   adapt(to, group.light_half_saturation),
                   adapt(to, group.light_limitation),
                   adapt(to, group.maximum_growth_rate),
                   adapt(to, group.mortality_rate),
                   adapt(to, group.quadratic_mortality_rate),
                   adapt(to, group.iron_ratio),
                   adapt(to, group.chlorophyll_ratio))

Base.summary(group::Photoautotroph) = string("Photoautotroph limited by $(limiting_nutrients(group))")

#####
##### Picoplankton and diazotrophs
#####

"""
    PhytoDiazotroph([FT = Float64;] kwargs...)
    PhytoDiazotroph(grid; picoplankton_sinking_speed = 0, diazotroph_sinking_speed = 0, open_bottom = true, kwargs...)

A plankton component for the `plankton` slot of a [`NutrientsPlanktonDetritus`](@ref) model which
resolves two photoautotrophs and no grazers: picoplankton (`P`), which grow on dissolved inorganic
nitrogen, and diazotrophs (`Diaz`), which instead fix dinitrogen. It is used by the [`WTSP`](@ref)
preset.

Both groups are light, phosphate, and iron limited, and picoplankton are additionally nitrogen
limited. Diazotrophs need no fixed nitrogen but their fixation is inhibited by dissolved inorganic
nitrogen and, being nitrogenase bearing, they carry a much larger iron quota than the picoplankton
(`diazotroph_iron_ratio` ≫ `picoplankton_iron_ratio`), which is what makes iron availability set
where they can grow. Grazing is not resolved, so mortality (which may be given a density dependent
part) transfers both groups to detritus.

A fraction `diazotroph_nitrogen_release_fraction` of gross fixation is released as dissolved
inorganic nitrogen rather than built into biomass, which is the pathway by which fixed nitrogen
fuels the picoplankton. This released nitrogen comes straight from N₂, so it carries no carbon,
phosphorus, or iron.

Conservation
============

Phosphorus, iron, carbon, and oxygen are conserved exactly. Nitrogen is conserved apart from the
gross nitrogen fixation source (see [`nitrogen_fixation`](@ref)), which is the only external supply
of any element in the model. Two routings make this exact:

- the diazotroph's *excess* iron (`diazotroph_iron_ratio - picoplankton_iron_ratio` per mmol N) is
  returned directly to the dissolved iron pool when diazotrophs die, so the single detritus pool
  stays at the model's one Fe:N ratio; and
- released fixed nitrogen enters only the nitrogen budget, so the carbon, phosphorus and iron taken
  up by diazotrophs scale with the fraction of fixation which becomes biomass.

Keyword Arguments
=================

- `picoplankton_nutrient_half_saturations`, `diazotroph_nutrient_half_saturations`: `NamedTuple`s of
  half-saturation constants keyed by limiting nutrient; their keys set which nutrients limit each
  group (mmol X / m³)
- `picoplankton_maximum_growth_rate`, `diazotroph_maximum_growth_rate`: gross growth rates (1 / s).
  Note that the diazotroph rate is *gross* of the nitrogen release, so their maximum rate of increase
  is `(1 - diazotroph_nitrogen_release_fraction) * diazotroph_maximum_growth_rate`
- `picoplankton_mortality_rate`, `diazotroph_mortality_rate`: linear mortality rates (1 / s)
- `picoplankton_quadratic_mortality_rate`, `diazotroph_quadratic_mortality_rate`: density dependent
  mortality rates, zero by default, which close the model where grazers are not resolved
  (1 / s / mmol N / m³)
- `picoplankton_iron_ratio`, `diazotroph_iron_ratio`: iron quotas (mol Fe / mol N)
- `picoplankton_chlorophyll_ratio`, `diazotroph_chlorophyll_ratio`: chlorophyll content
  (g Chl / mol N)
- `picoplankton_light_half_saturation`, `diazotroph_light_half_saturation`: light half-saturations
  (W / m²)
- `light_limitation`: the light-limitation functional form, [`ExponentialLightLimitation`](@ref) by
  default
- `diazotroph_nitrogen_release_fraction`: the fraction of gross fixation released as dissolved
  inorganic nitrogen rather than built into biomass
- `nitrogen_fixation_inhibition_half_saturation`: the dissolved inorganic nitrogen concentration at
  which fixation is halved (mmol N / m³)
- `carbon_ratio`, `phosphate_ratio`: the carbon and phosphorus content of all organic matter, which
  is at a fixed stoichiometry (mol X / mol N)
- `picoplankton_sinking_velocity`, `diazotroph_sinking_velocity`: sinking velocity fields, most
  easily set with the `grid` method's `picoplankton_sinking_speed` and `diazotroph_sinking_speed`,
  which are positive downwards, so a negative speed makes buoyant `Trichodesmium` rise (m / s)
"""
struct PhytoDiazotroph{PP, DZ, FT, PS, DS}
                                 picoplankton :: PP
                                   diazotroph :: DZ

         diazotroph_nitrogen_release_fraction :: FT
 nitrogen_fixation_inhibition_half_saturation :: FT

                                 carbon_ratio :: FT
                              phosphate_ratio :: FT

                picoplankton_sinking_velocity :: PS
                  diazotroph_sinking_velocity :: DS
end

function PhytoDiazotroph(FT = Float64;
                         picoplankton_nutrient_half_saturations = (nitrate = 0.1,     # mmol N / m³
                                                                   phosphate = 0.03,  # mmol P / m³
                                                                   iron = 8e-5),      # mmol Fe / m³
                         picoplankton_light_half_saturation = 2.17,                   # W / m²
                         picoplankton_maximum_growth_rate = 0.8 / day,                # 1 / s
                         picoplankton_mortality_rate = 0.06 / day,                    # 1 / s
                         picoplankton_quadratic_mortality_rate = 0,                   # 1 / s / mmol N / m³
                         picoplankton_iron_ratio = 3.3125e-5,                         # mol Fe / mol N (Fe:C = 5 μmol/mol)
                         picoplankton_chlorophyll_ratio = 1.137,                      # g Chl / mol N (C:Chl = 70 g/g)

                         diazotroph_nutrient_half_saturations = (phosphate = 0.05,    # mmol P / m³
                                                                 iron = 2e-3),        # mmol Fe / m³
                         diazotroph_light_half_saturation = 2.17,                     # W / m²
                         diazotroph_maximum_growth_rate = 0.2 / 0.7 / day,            # 1 / s (0.2/day net of release)
                         diazotroph_mortality_rate = 0.1 / day,                       # 1 / s
                         diazotroph_quadratic_mortality_rate = 0,                     # 1 / s / mmol N / m³
                         diazotroph_iron_ratio = 2.65e-4,                             # mol Fe / mol N (Fe:C = 40 μmol/mol)
                         diazotroph_chlorophyll_ratio = 1.137,                        # g Chl / mol N

                         light_limitation = ExponentialLightLimitation(),

                         diazotroph_nitrogen_release_fraction = 0.3,
                         nitrogen_fixation_inhibition_half_saturation = 0.3,          # mmol N / m³

                         carbon_ratio = 6.625,                                        # mol C / mol N
                         phosphate_ratio = 1/16,                                      # mol P / mol N

                         picoplankton_sinking_velocity = (u = ZeroField(), v = ZeroField(), w = ZeroField()),
                         diazotroph_sinking_velocity = (u = ZeroField(), v = ZeroField(), w = ZeroField()))

    picoplankton = Photoautotroph(FT;
                                  nutrient_half_saturations = picoplankton_nutrient_half_saturations,
                                  light_half_saturation = picoplankton_light_half_saturation,
                                  light_limitation,
                                  maximum_growth_rate = picoplankton_maximum_growth_rate,
                                  mortality_rate = picoplankton_mortality_rate,
                                  quadratic_mortality_rate = picoplankton_quadratic_mortality_rate,
                                  iron_ratio = picoplankton_iron_ratio,
                                  chlorophyll_ratio = picoplankton_chlorophyll_ratio)

    diazotroph = Photoautotroph(FT;
                                nutrient_half_saturations = diazotroph_nutrient_half_saturations,
                                light_half_saturation = diazotroph_light_half_saturation,
                                light_limitation,
                                maximum_growth_rate = diazotroph_maximum_growth_rate,
                                mortality_rate = diazotroph_mortality_rate,
                                quadratic_mortality_rate = diazotroph_quadratic_mortality_rate,
                                iron_ratio = diazotroph_iron_ratio,
                                chlorophyll_ratio = diazotroph_chlorophyll_ratio)

    return PhytoDiazotroph(picoplankton,
                           diazotroph,
                           convert(FT, diazotroph_nitrogen_release_fraction),
                           convert(FT, nitrogen_fixation_inhibition_half_saturation),
                           convert(FT, carbon_ratio),
                           convert(FT, phosphate_ratio),
                           picoplankton_sinking_velocity,
                           diazotroph_sinking_velocity)
end

function PhytoDiazotroph(grid::AbstractGrid{FT};
                         picoplankton_sinking_speed = zero(FT),
                         diazotroph_sinking_speed = zero(FT),
                         open_bottom = true,
                         kwargs...) where FT

    sinking_velocities = setup_velocity_fields((P = picoplankton_sinking_speed,
                                                Diaz = diazotroph_sinking_speed),
                                               grid, open_bottom; three_D = true)

    return PhytoDiazotroph(FT; picoplankton_sinking_velocity = sinking_velocities.P,
                               diazotroph_sinking_velocity = sinking_velocities.Diaz,
                               kwargs...)
end

const PhytoDiazotroph_NPD{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:PhytoDiazotroph}

required_biogeochemical_tracers(::PhytoDiazotroph) = (:P, :Diaz)
required_biogeochemical_auxiliary_fields(::PhytoDiazotroph) = (:PAR, )

biogeochemical_drift_velocity(bgc::PhytoDiazotroph_NPD, ::Val{:P}) =
    bgc.plankton.picoplankton_sinking_velocity

biogeochemical_drift_velocity(bgc::PhytoDiazotroph_NPD, ::Val{:Diaz}) =
    bgc.plankton.diazotroph_sinking_velocity

# the model's currency is nitrogen, and carbon and phosphorus are at a fixed stoichiometry
# everywhere. Iron is *not*, since the diazotroph's quota is larger, so `iron_ratio` is the ratio of
# everything else (picoplankton, detritus, and the waste routed through it)
@inline nitrogen_ratio(::PhytoDiazotroph, ::NPD{FT}) where FT = one(FT)
@inline carbon_ratio(plankton::PhytoDiazotroph, ::NPD) = plankton.carbon_ratio
@inline phosphate_ratio(plankton::PhytoDiazotroph, ::NPD) = plankton.phosphate_ratio
@inline iron_ratio(plankton::PhytoDiazotroph, ::NPD) = plankton.picoplankton.iron_ratio

"""The iron the diazotroph carries in excess of everything else, per mmol N (mol Fe / mol N)."""
@inline excess_diazotroph_iron_ratio(plankton::PhytoDiazotroph) =
    plankton.diazotroph.iron_ratio - plankton.picoplankton.iron_ratio

@inline chlorophyll(plankton::PhytoDiazotroph, model) =
    plankton.picoplankton.chlorophyll_ratio * model.tracers.P +
    plankton.diazotroph.chlorophyll_ratio * model.tracers.Diaz

#####
##### Growth and nitrogen fixation
#####

"""Dissolved inorganic nitrogen (mmol N / m³), which inhibits nitrogen fixation."""
@inline dissolved_inorganic_nitrogen(i, j, k, grid, ::Nutrients{<:SingleTracerNutrient}, ::NPD, fields) =
    @inbounds fields.N[i, j, k]

@inline dissolved_inorganic_nitrogen(i, j, k, grid, ::Nutrients{<:NitrateAmmonia}, ::NPD, fields) =
    @inbounds fields.NO₃[i, j, k] + fields.NH₄[i, j, k]

@inline dissolved_inorganic_nitrogen(i, j, k, grid, ::Nutrients, ::NPD{FT}, fields) where FT = zero(FT)

@inline picoplankton_growth(i, j, k, grid, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) =
    @inbounds specific_growth_rate(i, j, k, grid, plankton.picoplankton, bgc, fields, auxiliary_fields) * fields.P[i, j, k]

# gross fixation: light, phosphate, and iron limited, and inhibited by the fixed nitrogen which is
# cheaper for the diazotroph to use than fixing its own
@inline function nitrogen_fixation(i, j, k, grid, plankton::PhytoDiazotroph, bgc::NPD, fields, auxiliary_fields)
    kᴺ = plankton.nitrogen_fixation_inhibition_half_saturation

    @inbounds Diaz = fields.Diaz[i, j, k]

    DIN = dissolved_inorganic_nitrogen(i, j, k, grid, bgc.nutrients, bgc, fields)

    μ = specific_growth_rate(i, j, k, grid, plankton.diazotroph, bgc, fields, auxiliary_fields)

    return μ * kᴺ / (kᴺ + DIN) * Diaz
end

"""The fixed nitrogen which becomes diazotroph biomass (mmol N / m³ / s)."""
@inline function diazotroph_growth(i, j, k, grid, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields)
    γ = plankton.diazotroph_nitrogen_release_fraction

    return (1 - γ) * nitrogen_fixation(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
end

@inline picoplankton_mortality(i, j, k, grid, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) =
    @inbounds mortality(plankton.picoplankton, fields.P[i, j, k])

@inline diazotroph_mortality(i, j, k, grid, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) =
    @inbounds mortality(plankton.diazotroph, fields.Diaz[i, j, k])

@inline (bgc::PhytoDiazotroph_NPD)(i, j, k, grid, ::Val{:P}, clock, fields, auxiliary_fields) = (
    picoplankton_growth(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  - picoplankton_mortality(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
)

@inline (bgc::PhytoDiazotroph_NPD)(i, j, k, grid, ::Val{:Diaz}, clock, fields, auxiliary_fields) = (
    diazotroph_growth(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  - diazotroph_mortality(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
)

#####
##### Nutrient uptake
#####

# the nitrogen built into new biomass, which sets the carbon and phosphorus uptake. Only the
# picoplankton's share of this comes out of the dissolved inorganic nitrogen pool
@inline nutrient_uptake(i, j, k, grid, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) = (
    picoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
  + diazotroph_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
)

@inline nutrient_uptake(i, j, k, grid, ::Val{:N}, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) =
    picoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

# where nitrate and ammonia are resolved the picoplankton's uptake is split between them in the same
# proportion as they contribute to the nitrogen limitation
@inline function nitrogen_uptake_shares(i, j, k, grid, group::Photoautotroph, ::NPD{FT}, fields) where FT
    kNO₃ = nutrient_half_saturations(group, Val(:NO₃))
    kNH₄ = nutrient_half_saturations(group, Val(:NH₄))

    @inbounds begin
        NO₃ = fields.NO₃[i, j, k] / kNO₃
        NH₄ = fields.NH₄[i, j, k] / kNH₄
    end

    total = NO₃ + NH₄ + eps(zero(FT))

    return NO₃ / total, NH₄ / total
end

@inline function nutrient_uptake(i, j, k, grid, ::Val{:NO₃}, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields)
    nitrate_share, _ = nitrogen_uptake_shares(i, j, k, grid, plankton.picoplankton, bgc, fields)

    return nitrate_share * picoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
end

@inline function nutrient_uptake(i, j, k, grid, ::Val{:NH₄}, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields)
    _, ammonia_share = nitrogen_uptake_shares(i, j, k, grid, plankton.picoplankton, bgc, fields)

    return ammonia_share * picoplankton_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
end

# new diazotroph biomass draws down its own, larger, iron quota
@inline nutrient_uptake(i, j, k, grid, ::Val{:Fe}, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) = (
    iron_ratio(i, j, k, grid, plankton, bgc, fields)
  * nutrient_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
  + excess_diazotroph_iron_ratio(plankton)
  * diazotroph_growth(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
)

#####
##### Waste routing
#####

# grazing is not resolved, so all mortality goes to detritus and nothing is remineralised straight
# out of living biomass
@inline inorganic_waste(i, j, k, grid, ::PhytoDiazotroph, ::NPD{FT}, args...) where FT = zero(FT)
@inline dissolved_waste(i, j, k, grid, ::PhytoDiazotroph, ::NPD{FT}, args...) where FT = zero(FT)

@inline solid_waste(i, j, k, grid, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) = (
    picoplankton_mortality(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
  + diazotroph_mortality(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
)

# the released fraction of gross fixation goes straight back to the dissolved inorganic nitrogen
# pool (as ammonia where nitrate and ammonia are resolved). Having come from N₂ it carries no
# carbon, phosphorus, or iron, so it appears in the nitrogen waste alone
@inline inorganic_nitrogen_waste(i, j, k, grid, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) =
    plankton.diazotroph_nitrogen_release_fraction *
    nitrogen_fixation(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

# detritus carries the picoplankton's iron quota, so the diazotroph's excess iron is returned
# directly to the dissolved pool when it dies. Without this the two quotas could not both be
# represented by one detritus pool without losing (or creating) iron
@inline inorganic_iron_waste(i, j, k, grid, plankton::PhytoDiazotroph, bgc, fields, auxiliary_fields) =
    excess_diazotroph_iron_ratio(plankton) *
    diazotroph_mortality(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

#####
##### Admin
#####

Adapt.adapt_structure(to, plankton::PhytoDiazotroph) =
    PhytoDiazotroph(adapt(to, plankton.picoplankton),
                    adapt(to, plankton.diazotroph),
                    adapt(to, plankton.diazotroph_nitrogen_release_fraction),
                    adapt(to, plankton.nitrogen_fixation_inhibition_half_saturation),
                    adapt(to, plankton.carbon_ratio),
                    adapt(to, plankton.phosphate_ratio),
                    adapt(to, plankton.picoplankton_sinking_velocity),
                    adapt(to, plankton.diazotroph_sinking_velocity))

Base.summary(::PhytoDiazotroph) = "PhytoDiazotroph (:P, :Diaz)"

function Base.show(io::IO, plankton::PhytoDiazotroph)
    msg = summary(plankton) * "\n"

    msg *= "├── C:N:P = $(plankton.carbon_ratio):1:$(plankton.phosphate_ratio)\n"
    msg *= "├── Fe:N = $(plankton.picoplankton.iron_ratio) (picoplankton), $(plankton.diazotroph.iron_ratio) (diazotroph)\n"
    msg *= "├── Picoplankton limited by: $(limiting_nutrients(plankton.picoplankton))\n"
    msg *= "└── Diazotrophs limited by: $(limiting_nutrients(plankton.diazotroph)), and inhibited by dissolved inorganic nitrogen"

    print(io, msg)

    return nothing
end

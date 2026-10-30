using OceanBioME: ImplicitSinking, compute_implicit_sinking!

import OceanBioME: implicit_sinking_production

"""
    ImplicitCalciumCarbonate(grid;
                             replicates = 1,
                             particulate_fraction = 0.67,
                             dissolution_lengths = (soft_CaCO₃ = 500, hard_CaCO₃ = 40000),
                             production_partitioning = (soft_CaCO₃ = 0.98, hard_CaCO₃ = 0.02),
                             hard_organic_carbon_ratio = 0.01 * 100.09 / 12.01,
                             open_bottom = true,
                             store_flux = false)

The inorganic carbon component for the `inorganic_carbon` slot of a [`NutrientsPlanktonDetritus`](@ref) model
whose calcium carbonate sinks **implicitly** (see [`ImplicitSinking`](@ref)), alongside dissolved inorganic
carbon (`DIC`) and alkalinity (`Alk`).

Like [`CarbonateSystem`](@ref) there is no calcium carbonate tracer, but rather than dissolving with the
organic detritus the calcium carbonate has its own dissolution length scale: every step, the fraction
`particulate_fraction` of the calcium carbonate the plankton release (`particulate_calcium_carbonate_production`
plus `biological_calcium_carbonate_dissolution`, which for plankton with no standing biomass, such as
[`ImplicitProductivity`](@ref), is the rain ratio × primary production) is sunk down each column and dissolves
over `dissolution_lengths`, returning `DIC` and `2×` that to `Alk`. The rest is returned where it is released,
whatever the split of the organic waste. The default is MARBL's, where `f_graze_CaCO3_remin = 0.33` of the
calcium carbonate dissolves where it is grazed.

Like MARBL's hard particulate organic carbon, the sinking calcium carbonate carries organic matter with it:
`hard_organic_carbon_ratio` (mol C / mol CaCO₃, default MARBL's `parm_hPOC_CaCO3_ratio` of 0.01 g C / g CaCO₃)
times the calcium carbonate flux through every face (of all the classes), converted to the plankton's units with
its `carbon_ratio`. Its production is taken out of the particulate production of [`DissolvedParticulate`](@ref)
detritus (in proportion to its `particulate_waste_partitioning`) and it remineralises with the calcium
carbonate's dissolution like that detritus, so it needs no tracer. With an open bottom its floor flux (the ratio
times the calcium carbonate floor flux) is left for a sediment along with the detritus' floor flux. Other detritus
ignore it, so should be used with `hard_organic_carbon_ratio = 0`, which turns it off.

The production is split between the classes of `production_partitioning`, which sink independently with the
`dissolution_lengths` of the same names (m, numbers or functions `(i, j, k, grid, clock, fields)` as for
`ImplicitSinking`). The defaults follow MARBL's base `parm_CaCO3_diss = 500 m` with the fraction
`parm_CaCO3_gamma = 0.02` of the production in a hard subclass which dissolves over 40 000 m, but not MARBL's
lengthening with depth and at low O₂.

With `open_bottom = true` the calcium carbonate reaching the floor leaves the domain, and its flux (the sum of
the `sinking.floor_flux` of the classes) is left for a sediment; otherwise it dissolves in the floor cell.
`store_flux = true` also stores the flux through each cell's bottom face (`sinking.flux`).

The production does not depend on the carbonate tracers, so the sinking is solved once and serves every
replicate: passing `replicates > 1` adds `replicates` copies of the carbonate tracers (`DIC1`, `Alk1`, `DIC2`, …)
which evolve with the same tendency as the base `DIC`/`Alk`, and differ only through what a sediment does with
the floor flux.
"""
struct ImplicitCalciumCarbonate{N, FT, PP, SK} <: AbstractInorganicCarbon{N}
         particulate_fraction :: FT # –
      production_partitioning :: PP # –
    hard_organic_carbon_ratio :: FT # mol C / mol CaCO₃
                      sinking :: SK

    ImplicitCalciumCarbonate{N}(particulate_fraction::FT, production_partitioning::PP,
                                hard_organic_carbon_ratio::FT, sinking::SK) where {N, FT, PP, SK} =
        new{N, FT, PP, SK}(particulate_fraction, production_partitioning, hard_organic_carbon_ratio, sinking)
end

function ImplicitCalciumCarbonate(grid::AbstractGrid{FT};
                                  replicates = 1,
                                  particulate_fraction = 0.67,
                                  dissolution_lengths = (soft_CaCO₃ = 500, hard_CaCO₃ = 40000),
                                  production_partitioning = (soft_CaCO₃ = 0.98, hard_CaCO₃ = 0.02),
                                  hard_organic_carbon_ratio = 0.01 * 100.09 / 12.01, # mol C / mol CaCO₃
                                  open_bottom = true,
                                  store_flux = false) where FT

    0 ≤ particulate_fraction ≤ 1 ||
        throw(ArgumentError("the `particulate_fraction` is $particulate_fraction, which is not between 0 and 1"))

    hard_organic_carbon_ratio ≥ 0 ||
        throw(ArgumentError("the `hard_organic_carbon_ratio` is $hard_organic_carbon_ratio, which is negative"))

    keys(production_partitioning) == keys(dissolution_lengths) ||
        throw(ArgumentError("the `production_partitioning` is for $(keys(production_partitioning)) but the " *
                            "`dissolution_lengths` are for $(keys(dissolution_lengths))"))

    sum(production_partitioning) ≈ 1 ||
        throw(ArgumentError("the `production_partitioning` sums to $(sum(production_partitioning)), not 1, " *
                            "so calcium carbonate would not be conserved"))

    production_partitioning = NamedTuple{keys(production_partitioning)}(convert.(FT, values(production_partitioning)))

    sinking = ImplicitSinking(grid, dissolution_lengths; open_bottom, store_flux)

    return ImplicitCalciumCarbonate{replicates}(convert(FT, particulate_fraction),
                                                production_partitioning,
                                                convert(FT, hard_organic_carbon_ratio),
                                                sinking)
end

const NPD_ICC{FT} = NPD{FT, <:Any, <:Any, <:Any, <:ImplicitCalciumCarbonate}

required_biogeochemical_tracers(::ImplicitCalciumCarbonate{1}) = (:DIC, :Alk)
required_biogeochemical_tracers(::ImplicitCalciumCarbonate{N}) where N =
    (map(n->Symbol(:DIC, n), 1:N)..., map(n->Symbol(:Alk, n), 1:N)...)

required_biogeochemical_auxiliary_fields(::ImplicitCalciumCarbonate) = tuple()

# the calcium carbonate the plankton release, which for plankton with no standing biomass (e.g. `ImplicitProductivity`)
# is all that they form
@inline released_calcium_carbonate(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
    particulate_calcium_carbonate_production(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) +
    biological_calcium_carbonate_dissolution(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

# `particulate_fraction` of it sinks and the rest dissolves where it is released
@inline sinking_calcium_carbonate_production(i, j, k, grid, ic::ImplicitCalciumCarbonate, bgc, fields, auxiliary_fields) =
    ic.particulate_fraction * released_calcium_carbonate(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)

@inline implicit_sinking_production(i, j, k, grid, ic::ImplicitCalciumCarbonate, bgc, fields, auxiliary_fields, ::Val{name}) where name =
    getproperty(ic.production_partitioning, name) *
    sinking_calcium_carbonate_production(i, j, k, grid, ic, bgc, fields, auxiliary_fields)

# the dissolution of all the classes, which includes what is formed in the cell and dissolves there
@inline sinking_calcium_carbonate_dissolution(i, j, k, remineralisation::Tuple{Any}) =
    @inbounds first(remineralisation)[i, j, k]

@inline sinking_calcium_carbonate_dissolution(i, j, k, remineralisation::Tuple) =
    @inbounds first(remineralisation)[i, j, k] + sinking_calcium_carbonate_dissolution(i, j, k, Base.tail(remineralisation))

@inline net_calcium_carbonate_production(i, j, k, grid, bgc::NPD_ICC, fields, auxiliary_fields) =
    biological_calcium_carbonate_precipitation(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) -
    (1 - bgc.inorganic_carbon.particulate_fraction) * released_calcium_carbonate(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) -
    sinking_calcium_carbonate_dissolution(i, j, k, values(bgc.inorganic_carbon.sinking.remineralisation))

# the hard organic matter carried per calcium carbonate, in the plankton's units (e.g. mol P / mol CaCO₃), so that
# its flux through each face is this times the calcium carbonate flux there
@inline hard_organic_ratio(i, j, k, grid, ic::ImplicitCalciumCarbonate, bgc, fields) =
    ic.hard_organic_carbon_ratio / carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields)

@inline hard_organic_production(i, j, k, grid, ic::ImplicitCalciumCarbonate, bgc::NPD, fields, auxiliary_fields) =
    hard_organic_ratio(i, j, k, grid, ic, bgc, fields) *
    sinking_calcium_carbonate_production(i, j, k, grid, ic, bgc, fields, auxiliary_fields)

@inline hard_organic_remineralisation(i, j, k, grid, ic::ImplicitCalciumCarbonate, bgc::NPD, fields, auxiliary_fields) =
    hard_organic_ratio(i, j, k, grid, ic, bgc, fields) *
    sinking_calcium_carbonate_dissolution(i, j, k, values(ic.sinking.remineralisation))

update_biogeochemical_state!(model, ic::ImplicitCalciumCarbonate, npd::NPD) =
    compute_implicit_sinking!(model, ic.sinking, ic, npd)

# With `replicates > 1` the tracers are `DIC1`, `Alk1`, `DIC2`, … rather than `DIC` and `Alk`, and
# each replicate evolves with the tendency of the base tracer (see `component_tendency` for `CarbonateSystem`).
@inline @generated function component_tendency(i, j, k, grid, ::ImplicitCalciumCarbonate{N}, ::Val{name},
                                               bgc::NPD_ICC{FT}, clock, fields, auxiliary_fields) where {N, name, FT}
    N > 1 || return :(zero($FT))

    for n in 1:N
        name == Symbol(:DIC, n) && return :($(Expr(:meta, :inline)); bgc(i, j, k, grid, Val(:DIC), clock, fields, auxiliary_fields))
        name == Symbol(:Alk, n) && return :($(Expr(:meta, :inline)); bgc(i, j, k, grid, Val(:Alk), clock, fields, auxiliary_fields))
    end

    return :(zero($FT))
end

Adapt.adapt_structure(to, ic::ImplicitCalciumCarbonate{N}) where N =
    ImplicitCalciumCarbonate{N}(adapt(to, ic.particulate_fraction),
                                adapt(to, ic.production_partitioning),
                                adapt(to, ic.hard_organic_carbon_ratio),
                                adapt(to, ic.sinking))

Base.summary(ic::ImplicitCalciumCarbonate{1}) =
    string("ImplicitCalciumCarbonate $(required_biogeochemical_tracers(ic))")

Base.summary(ic::ImplicitCalciumCarbonate{N}) where N =
    string("ImplicitCalciumCarbonate{realisations = $N} $(required_biogeochemical_tracers(ic))")

function Base.show(io::IO, ic::ImplicitCalciumCarbonate{N}) where N
    msg = summary(ic) * "\n"

    if N > 1
        msg *= "├── Realisations: $N\n"
    end

    msg *= "├── Particulate fraction: $(ic.particulate_fraction)\n"
    msg *= "├── Production partitioning: $(ic.production_partitioning)\n"
    msg *= "├── Hard organic carbon ratio: $(ic.hard_organic_carbon_ratio) mol C / mol CaCO₃\n"
    msg *= "└── Sinking: $(summary(ic.sinking))"

    print(io, msg)

    return nothing
end

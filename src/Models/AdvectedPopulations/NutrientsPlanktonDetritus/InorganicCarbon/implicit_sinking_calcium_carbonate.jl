using OceanBioME: ImplicitSinking, compute_implicit_sinking!

import OceanBioME: implicit_sinking_production

"""
    ImplicitCalciumCarbonate(grid;
                             replicates = 1,
                             dissolution_lengths = (soft_CaCO₃ = 500, hard_CaCO₃ = 40000),
                             production_partitioning = (soft_CaCO₃ = 0.98, hard_CaCO₃ = 0.02),
                             open_bottom = true,
                             store_flux = false)

The inorganic carbon component for the `inorganic_carbon` slot of a [`NutrientsPlanktonDetritus`](@ref) model
whose calcium carbonate sinks **implicitly** (see [`ImplicitSinking`](@ref)), alongside dissolved inorganic
carbon (`DIC`) and alkalinity (`Alk`).

Like [`CarbonateSystem`](@ref) there is no calcium carbonate tracer, but rather than dissolving with the
organic detritus the calcium carbonate has its own dissolution length scale: every step, the calcium carbonate
formed with the plankton's solid waste (`particulate_calcium_carbonate_production`) is sunk down each column and
dissolves over `dissolution_lengths`, returning `DIC` and `2×` that to `Alk`. The calcium carbonate formed with
the dissolved and inorganic waste (`biological_calcium_carbonate_dissolution`) is returned where it forms.

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
struct ImplicitCalciumCarbonate{N, PP, SK} <: AbstractInorganicCarbon{N}
    production_partitioning :: PP # –
                    sinking :: SK

    ImplicitCalciumCarbonate{N}(production_partitioning::PP, sinking::SK) where {N, PP, SK} =
        new{N, PP, SK}(production_partitioning, sinking)
end

function ImplicitCalciumCarbonate(grid::AbstractGrid{FT};
                                  replicates = 1,
                                  dissolution_lengths = (soft_CaCO₃ = 500, hard_CaCO₃ = 40000),
                                  production_partitioning = (soft_CaCO₃ = 0.98, hard_CaCO₃ = 0.02),
                                  open_bottom = true,
                                  store_flux = false) where FT

    keys(production_partitioning) == keys(dissolution_lengths) ||
        throw(ArgumentError("the `production_partitioning` is for $(keys(production_partitioning)) but the " *
                            "`dissolution_lengths` are for $(keys(dissolution_lengths))"))

    sum(production_partitioning) ≈ 1 ||
        throw(ArgumentError("the `production_partitioning` sums to $(sum(production_partitioning)), not 1, " *
                            "so calcium carbonate would not be conserved"))

    production_partitioning = NamedTuple{keys(production_partitioning)}(convert.(FT, values(production_partitioning)))

    sinking = ImplicitSinking(grid, dissolution_lengths; open_bottom, store_flux)

    return ImplicitCalciumCarbonate{replicates}(production_partitioning, sinking)
end

const NPD_ICC{FT} = NPD{FT, <:Any, <:Any, <:Any, <:ImplicitCalciumCarbonate}

required_biogeochemical_tracers(::ImplicitCalciumCarbonate{1}) = (:DIC, :Alk)
required_biogeochemical_tracers(::ImplicitCalciumCarbonate{N}) where N =
    (map(n->Symbol(:DIC, n), 1:N)..., map(n->Symbol(:Alk, n), 1:N)...)

required_biogeochemical_auxiliary_fields(::ImplicitCalciumCarbonate) = tuple()

@inline implicit_sinking_production(i, j, k, grid, ic::ImplicitCalciumCarbonate, bgc, fields, auxiliary_fields, ::Val{name}) where name =
    getproperty(ic.production_partitioning, name) *
    particulate_calcium_carbonate_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)

# the dissolution of all the classes, which includes what is formed in the cell and dissolves there
@inline sinking_calcium_carbonate_dissolution(i, j, k, remineralisation::Tuple{Any}) =
    @inbounds first(remineralisation)[i, j, k]

@inline sinking_calcium_carbonate_dissolution(i, j, k, remineralisation::Tuple) =
    @inbounds first(remineralisation)[i, j, k] + sinking_calcium_carbonate_dissolution(i, j, k, Base.tail(remineralisation))

@inline net_calcium_carbonate_production(i, j, k, grid, bgc::NPD_ICC, fields, auxiliary_fields) =
    biological_calcium_carbonate_precipitation(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) -
    biological_calcium_carbonate_dissolution(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) -
    sinking_calcium_carbonate_dissolution(i, j, k, values(bgc.inorganic_carbon.sinking.remineralisation))

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
    ImplicitCalciumCarbonate{N}(adapt(to, ic.production_partitioning),
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

    msg *= "├── Production partitioning: $(ic.production_partitioning)\n"
    msg *= "└── Sinking: $(summary(ic.sinking))"

    print(io, msg)

    return nothing
end

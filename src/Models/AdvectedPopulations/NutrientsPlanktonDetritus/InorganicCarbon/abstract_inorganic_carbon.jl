# An inorganic carbon component with a replicate for each of the labels `L`, each of which has its own
# `tracers` (e.g. `(:DIC, :Alk)`) named by `replicate_name`
abstract type AbstractInorganicCarbon{L, tracers} end

replicate_labels(replicates::Integer) = replicates == 1 ? (nothing, ) : Tuple(1:replicates)
replicate_labels(label::Symbol) = (label, )

function replicate_labels(labels)
    allunique(labels) || throw(ArgumentError("The replicate labels $labels are not unique"))
    return Tuple(labels)
end

replicate_name(name, ::Nothing) = name
replicate_name(name, n::Integer) = Symbol(name, n)
replicate_name(name, label::Symbol) = Symbol(name, :_, label)

# e.g. `fields.CaCO₃_control` for `replicate_field(fields, Val(:CaCO₃), Val(:control))`, resolved at compile time
@inline @generated replicate_field(fields, ::Val{name}, ::Val{label}) where {name, label} =
    :(fields.$(replicate_name(name, label)))

required_biogeochemical_tracers(::AbstractInorganicCarbon{L, tracers}) where {L, tracers} =
    Tuple(replicate_name(tracer, label) for tracer in tracers for label in L)

carbonate_fields(::AbstractInorganicCarbon{L, tracers}) where {L, tracers} = tracers

"""
    carbonate_replicates(bgc)

Return the label and the field names of each replicate of the inorganic carbon component of `bgc`,
which may be a [`NutrientsPlanktonDetritus`](@ref) model, a `Biogeochemistry` wrapping one, or the
component itself. The single unlabelled replicate (e.g. `CarbonateSystem()`) has the label `nothing`
and numbered replicates (e.g. `CarbonateSystem(2)`) have the labels `1`, `2`, ….

```jldoctest
julia> using OceanBioME

julia> carbonate_replicates(CarbonateSystem((:control, :release)))
((label = :control, DIC = :DIC_control, Alk = :Alk_control), (label = :release, DIC = :DIC_release, Alk = :Alk_release))
```
"""
carbonate_replicates(ic::AbstractInorganicCarbon{L}) where L =
    map(label -> (; label, NamedTuple{carbonate_fields(ic)}(replicate_name.(carbonate_fields(ic), label))...), L)

carbonate_replicates(bgc::NutrientsPlanktonDetritus) = carbonate_replicates(bgc.inorganic_carbon)
carbonate_replicates(bgc::CompleteBiogeochemistry) = carbonate_replicates(bgc.underlying_biogeochemistry)
carbonate_replicates(::Nothing) = ()

# e.g. a sediment coupled to `DIC` exchanges it with the `DIC` of every replicate
tracer_replicates(bgc::NutrientsPlanktonDetritus, name) = tracer_replicates(bgc.inorganic_carbon, name)
tracer_replicates(::AbstractInorganicCarbon{L, tracers}, name) where {L, tracers} =
    name in tracers ? replicate_name.(name, L) : (name, )

Base.summary(ic::AbstractInorganicCarbon{L}) where L =
    string(nameof(typeof(ic)), length(L) > 1 ? "{realisations = $(length(L))} " : " ", required_biogeochemical_tracers(ic))

@inline net_biological_dic_uptake(i, j, k, grid, bgc, fields, auxiliary_fields) = (
  - primary_production(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_carbon_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_carbon_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
)

@inline function net_biological_alkalinity_uptake(i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    PN = phosphate_ratio(i, j, k, grid, bgc.plankton, bgc, fields) / nitrogen_ratio(i, j, k, grid, bgc.plankton, bgc, fields)

    return (
      - nitrate_primary_production(i, j, k, grid, bgc.plankton, bgc, clock, fields, auxiliary_fields) * (1 + PN)
      + ammonia_primary_production(i, j, k, grid, bgc.plankton, bgc, clock, fields, auxiliary_fields) * (1 - PN)
    )
end

@inline carbonate_tendency(::Val{:DIC}, i, j, k, grid, bgc, clock, fields, auxiliary_fields, replicate, label) =
    net_biological_dic_uptake(i, j, k, grid, bgc, fields, auxiliary_fields) -
    net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields, replicate, label)

@inline carbonate_tendency(::Val{:Alk}, i, j, k, grid, bgc, clock, fields, auxiliary_fields, replicate, label) =
    net_biological_alkalinity_uptake(i, j, k, grid, bgc, clock, fields, auxiliary_fields) -
    2 * net_calcium_carbonate_production(i, j, k, grid, bgc, fields, auxiliary_fields, replicate, label)

# The tracers of the replicates are only named when the component is constructed (from the labels `L`
# and `tracers`), so `name` is matched against them here at compile time and gets the tendency of its
# kind of tracer in its replicate (index and label). A name no replicate has, such as a passive `DIC`
# beside `CarbonateSystem(2)`, gets nothing (see `component_tendency` in nutrients_plankton_detritus.jl).
@inline @generated function component_tendency(i, j, k, grid, ::AbstractInorganicCarbon{L, tracers}, ::Val{name},
                                               bgc::NPD{FT}, clock, fields, auxiliary_fields) where {L, tracers, name, FT}
    for (n, label) in enumerate(L), tracer in tracers
        name == replicate_name(tracer, label) && return quote
            $(Expr(:meta, :inline))
            carbonate_tendency($(Val(tracer)), i, j, k, grid, bgc, clock, fields, auxiliary_fields, $(Val(n)), $(Val(label)))
        end
    end

    return :(zero($FT))
end

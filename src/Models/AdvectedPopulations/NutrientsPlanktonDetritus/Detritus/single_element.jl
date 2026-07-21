"""
    DissolvedParticulate(grid, dissolved_names = :DOM, particulate_names = (:sPOM, :bPOM); kwargs...)

A single-element detritus component for the `detritus` slot of a [`NutrientsPlanktonDetritus`](@ref)
model that splits organic matter into any number of dissolved classes (which do not sink) and
particulate classes (which sink). The tracers added are the `dissolved_names` and `particulate_names`
(e.g. the default `DOM`, `sPOM`, `bPOM`). Plankton dissolved/solid waste is partitioned between the
classes, particulate matter sinks and remineralises, and a fraction of particulate remineralisation
may be routed through the dissolved pool.

Keyword Arguments
=================

- `grid`: (required) the geometry, needed to configure the sinking-speed fields
- `dissolved_names`, `particulate_names`: the tracer name(s) of each dissolved / particulate class
- `dissolved_remineralisation_rate`, `particulate_remineralisation_rate`: per-class remineralisation
  rates (1/s)
- `dissolved_waste_partitioning`, `particulate_waste_partitioning`: how dissolved / solid plankton
  waste is split across the classes (each sums to 1)
- `dissolved_fraction_of_remineralisation`: the fraction of each particulate class's remineralisation
  that passes through the dissolved pool rather than directly to inorganic nutrients
- `sinking_speeds`: the downward sinking speed of each particulate class (m/s)
- `open_bottom`: whether particulate detritus can sink out of the bottom of the domain
"""
struct DissolvedParticulate{N, M, DN, PN, FN, FM, SV}
         dissolved_remineralisation_rate :: FN
       particulate_remineralisation_rate :: FM
        
            dissolved_waste_partitioning :: FN
          particulate_waste_partitioning :: FM

  dissolved_fraction_of_remineralisation :: FM

                      sinking_velocities :: SV
end

function DissolvedParticulate(FT = Float64;
                              dissolved_remineralisation_rate,
                              particulate_remineralisation_rate,
                              dissolved_waste_partitioning,
                              particulate_waste_partitioning,
                              dissolved_fraction_of_remineralisation,
                              sinking_velocities::SV,
                              dissolved_names,
                              particulate_names) where SV

    dissolved_names = possibly_tuple_or_symbol(dissolved_names)
    particulate_names = possibly_tuple_or_symbol(particulate_names)

    dissolved_remineralisation_rate = convert.(FT, dissolved_remineralisation_rate)
    particulate_remineralisation_rate = convert.(FT, particulate_remineralisation_rate)
    dissolved_waste_partitioning = convert.(FT, dissolved_waste_partitioning)
    particulate_waste_partitioning = convert.(FT, particulate_waste_partitioning)
    dissolved_fraction_of_remineralisation = convert.(FT, dissolved_fraction_of_remineralisation)

    DN = typeof(dissolved_names)
    PN = typeof(particulate_names)
    FN = typeof(dissolved_remineralisation_rate)
    FM = typeof(particulate_remineralisation_rate)

    return DissolvedParticulate{length(dissolved_names), 
                                length(particulate_names),
                                dissolved_names,
                                particulate_names,
                                FN, FM, SV}(dissolved_remineralisation_rate, 
                                            particulate_remineralisation_rate,
                                            dissolved_waste_partitioning,
                                            particulate_waste_partitioning,
                                            dissolved_fraction_of_remineralisation,
                                            sinking_velocities)
end

required_biogeochemical_tracers(dp::DissolvedParticulate{N, M, DN, PN}) where {N, M, DN, PN} = 
    (DN..., PN...)

required_biogeochemical_auxiliary_fields(::DissolvedParticulate) = tuple()

const NPD_DP{FT, N, M} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:DissolvedParticulate{N, M}}

function DissolvedParticulate(grid::AbstractGrid{FT}, dissolved_names = :DOM, particulate_names = (:sPOM, :bPOM);
                              dissolved_remineralisation_rate = repeat_property(dissolved_names, 3.86e-7), # 1/s
                              particulate_remineralisation_rate =  repeat_property(particulate_names, 5.88e-7),
                              dissolved_waste_partitioning =  default_partitioning(dissolved_names),
                              particulate_waste_partitioning = default_partitioning(particulate_names),
                              dissolved_fraction_of_remineralisation = repeat_property(particulate_names, one(FT)),
                              sinking_speeds = default_sinking_speeds(particulate_names),
                              open_bottom = true) where FT

    sinking_velocities = setup_velocity_fields(NamedTuple{possibly_tuple_or_symbol(particulate_names)}(sinking_speeds), grid, open_bottom; three_D = true)

    manifest_multi_class_dissolved_particulate(dissolved_names, particulate_names)

    return DissolvedParticulate(FT;
                                dissolved_remineralisation_rate,
                                particulate_remineralisation_rate,
                                dissolved_waste_partitioning,
                                particulate_waste_partitioning,
                                dissolved_fraction_of_remineralisation,
                                sinking_velocities,
                                dissolved_names,
                                particulate_names)
end

possibly_tuple_or_symbol(names) = names
possibly_tuple_or_symbol(names::Symbol) = tuple(names)
repeat_property(names, value) = tuple(repeat([value], length(names))...)
repeat_property(::Symbol, value) = value

default_sinking_speeds(names, FT=Float64) = tuple(repeat([convert(FT, 10/day)], length(names))...)
default_sinking_speeds(::NTuple{2}, FT=Float64) = (convert(FT, 3/day), convert(FT, 200/day))
default_sinking_speeds(::Symbol, FT=Float64) = convert(FT, 10/day)

default_partitioning(names, FT=Float64) = tuple(repeat([convert(FT, 1/length(names))], length(names))...)
default_partitioning(::Symbol, FT=Float64) = one(FT)

const _manifested_dissolved_particulate = Set{Tuple}()

function manifest_multi_class_dissolved_particulate(dissolved_names, particulate_names)
    dissolved_names = possibly_tuple_or_symbol(dissolved_names)
    particulate_names = possibly_tuple_or_symbol(particulate_names)

    key = (dissolved_names, particulate_names)
    key in _manifested_dissolved_particulate && return nothing
    push!(_manifested_dissolved_particulate, key)

    for (n, name) in enumerate(dissolved_names)
        @eval begin
            @inline (bgc::NPD_DP)(i, j, k, grid, val_name::Val{$(QuoteNode(name))}, clock, fields, auxiliary_fields) = @inbounds (
                dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) * bgc.detritus.dissolved_waste_partitioning[$n]
              + dissolved_remineralisation(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields) * bgc.detritus.dissolved_waste_partitioning[$n]
              - grazing(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
              - bgc.detritus.dissolved_remineralisation_rate[$n] * fields[$(QuoteNode(name))][i, j, k]
            )
        end
    end

    for (m, name) in enumerate(particulate_names)
        @eval begin
            @inline (bgc::NPD_DP)(i, j, k, grid, val_name::Val{$(QuoteNode(name))}, clock, fields, auxiliary_fields) = @inbounds (
                solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields) * bgc.detritus.particulate_waste_partitioning[$m]
              - grazing(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
              - bgc.detritus.particulate_remineralisation_rate[$m] * fields[$(QuoteNode(name))][i, j, k]
            )

            @inline biogeochemical_drift_velocity(bgc::NPD_DP, ::Val{$(QuoteNode(name))}) =
                bgc.detritus.sinking_velocities[$m]
        end
    end

    return nothing
end

@inline @generated function dissolved_remineralisation(i, j, k, grid, detritus::DissolvedParticulate{N, M, DN, PN}, bgc::NPD_DP{FT}, fields, auxiliary_fields) where {N, M, DN, PN, FT}
    combined = Expr(:block)
    push!(combined.args, :(total = zero($FT)))
    for (m, name) in enumerate(PN)
        push!(combined.args, :(total += fields[$(QuoteNode(name))][i, j, k] *
                        detritus.particulate_remineralisation_rate[$m] *
                        detritus.dissolved_fraction_of_remineralisation[$m]))
    end
    push!(combined.args, :(return total))
    return combined
end

@inline @generated function inorganic_waste(i, j, k, grid, detritus::DissolvedParticulate{N, M, DN, PN}, bgc::NPD_DP{FT}, fields, auxiliary_fields) where {N, M, DN, PN, FT}
    combined = Expr(:block)
    push!(combined.args, :(total = zero($FT)))
    for (m, name) in enumerate(PN)
        push!(combined.args, :(total += fields[$(QuoteNode(name))][i, j, k] *
                        detritus.particulate_remineralisation_rate[$m] *
                        (one($FT) - detritus.dissolved_fraction_of_remineralisation[$m])))
    end
    for (n, name) in enumerate(DN)
        push!(combined.args, :(total += fields[$(QuoteNode(name))][i, j, k] *
                        detritus.dissolved_remineralisation_rate[$n]))
    end
    push!(combined.args, :(return total))
    return combined
end

@inline @generated function calcite_dissolution(i, j, k, grid, detritus::DissolvedParticulate{N, M, DN, PN}, bgc::NPD_DP{FT}, fields, auxiliary_fields) where {N, M, DN, PN, FT}
    combined = Expr(:block)
    push!(combined.args, :(total = zero($FT)))
    for (m, name) in enumerate(PN)
        push!(combined.args, :(total += fields[$(QuoteNode(name))][i, j, k] *
                        detritus.particulate_remineralisation_rate[$m]))
    end
    push!(combined.args, :(total += dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)))
    push!(combined.args, :(return total * carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields) *
                                         calcite_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)))
    return combined
end

# admin
function Adapt.adapt_structure(to, detritus::DissolvedParticulate{N, M, DN, PN}) where {N, M, DN, PN}
    dissolved_remineralisation_rate = adapt(to, detritus.dissolved_remineralisation_rate)
    particulate_remineralisation_rate = adapt(to, detritus.particulate_remineralisation_rate)
    dissolved_waste_partitioning = adapt(to, detritus.dissolved_waste_partitioning)
    particulate_waste_partitioning = adapt(to, detritus.particulate_waste_partitioning)
    dissolved_fraction_of_remineralisation = adapt(to, detritus.dissolved_fraction_of_remineralisation)
    sinking_velocities = adapt(to, detritus.sinking_velocities)

    FN = typeof(dissolved_remineralisation_rate)
    FM = typeof(particulate_remineralisation_rate)
    SV = typeof(sinking_velocities)
    
    return DissolvedParticulate{N, M, DN, PN, FN, FM, SV}(
        dissolved_remineralisation_rate,
        particulate_remineralisation_rate,
        dissolved_waste_partitioning,
        particulate_waste_partitioning,
        dissolved_fraction_of_remineralisation,
        sinking_velocities # I don't know if we need this
    )
end

Base.summary(dp::DissolvedParticulate{N, M}) where {N, M} = 
    string("DissolvedParticulate{dissolved=$N, particulate=$M} $(required_biogeochemical_tracers(dp))")

Base.summary(dp::DissolvedParticulate{1, 1}) = 
    string("DissolvedParticulate $(required_biogeochemical_tracers(dp))")

Base.summary(dp::DissolvedParticulate{N, 1}) where N = 
    string("DissolvedParticulate{dissolved=$N} $(required_biogeochemical_tracers(dp))")

Base.summary(dp::DissolvedParticulate{1, M}) where M = 
    string("DissolvedParticulate{particulate=$M} $(required_biogeochemical_tracers(dp))")

function Base.show(io::IO, dp::DissolvedParticulate{N, M, ND, NP}) where {N, M, ND, NP}
    msg = summary(dp) * "\n"

    if N>1
        msg *= "├── Dissolved waste partitioning : $(dp.dissolved_waste_partitioning)\n"
    end

    if M>1
        msg *= "├── Particulate waste partitioning : $(dp.particulate_waste_partitioning)\n"
    end

    msg *= "└── Particle sinking speeds\n"

    for m in 1:M-1
        msg *= "  ├── $(NP[m]) : " * summary(dp.sinking_velocities[m].w) * "\n"
    end

    msg *= "  └── $(NP[end]) : " * summary(dp.sinking_velocities[end].w)

    print(io, msg)

    return nothing
end

using Oceananigans: fields
using KernelAbstractions: @kernel, @index
using Oceananigans.Architectures: architecture
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.DistributedComputations: synchronize_communication!
using Oceananigans.ImmersedBoundaries: immersed_cell
using Oceananigans.Utils: launch!

import Base: summary, show
import Oceananigans.Biogeochemistry: update_tendencies!, update_biogeochemical_state!

"""
    ZeroNegativeTracers(; exclude = ())

Construct a modifier that zeroes any negative tracers excluding those listed in `exclude`.

!!! danger "Tracer conservation"
    This method is _not_ recommended as a way to preserve positivity of tracers since
    it does not conserve the total tracer.
"""
@kwdef struct ZeroNegativeTracers{E}
    exclude :: E = ()
end

function update_biogeochemical_state!(model, zero::ZeroNegativeTracers)
    for (tracer_name, tracer) in pairs(model.tracers)
        if !(tracer_name in zero.exclude)
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end

#####
##### Infastructure to rescale negative values
#####

struct ScaleNegativeTracers{FA, SA, SH, OT, FV}
             tracers :: FA
        scalefactors :: SA
              shared :: SH
    first_own_tracer :: OT
  invalid_fill_value :: FV
end

"""
    ScaleNegativeTracers(tracers; scalefactors = ones(length(tracers)), invalid_fill_value = NaN)
    ScaleNegativeTracers(groups::NamedTuple; invalid_fill_value = NaN)

Constructs a modifier to scale `tracers` so that none are negative, while conserving their total
weighted by `scalefactors`. Use like:
```julia
modifier = ScaleNegativeTracers((:P, :Z, :N))
biogeochemistry = Biogeochemistry(...; modifier)
```
This method is better, though still imperfect, method to prevent numerical errors that lead to
negative tracer values compared to [`ZeroNegativeTracers`](@ref). Please see [discussion in
github](https://github.com/OceanBioME/OceanBioME.jl/discussions/48).

In each cell where a tracer is negative, the negative tracers are set to zero and the positive ones are
multiplied by `t / p`, where `t` is the weighted total and `p` is the weighted total of the positive tracers.

Several totals can be conserved at once, for example the nitrogen and the carbon in the same plankton, by
passing a named tuple of groups, each a tuple of tracer names or a named tuple of scale factors:
```julia
modifier = ScaleNegativeTracers((nitrogen = (N = 1, P = 1, Z = 1),
                                 carbon = (DIC = 1, P = 6.56, Z = 6.56)))
```
A tracer in more than one group (`P` and `Z` here) is scaled by the smallest of its groups' factors `t / p`,
and the tracers that are only in one group (`N` and `DIC`) are then scaled to keep that group's total, so
every group is conserved (if none of them is positive, the first one gets the difference). Groups which
share tracers must therefore be given to the same `ScaleNegativeTracers`, and each group needs at least one
tracer of its own. Scale factors can not be negative, and a tracer with a zero scale factor is not part of
the group.

After scaling, the halos of the tracers are filled again so that the fluxes through the boundaries use the
scaled values.

Future plans include implement a positivity-preserving timestepping scheme as the ideal alternative.

~~If `warn` is true then scaling will raise a warning.~~

`invalid_fill_value` specifies the value to set the total cell content to if the total is less than 0
(meaning that total tracer conservation can not be enforced). If the value is set to anything other than
`NaN` this scheme no longer conserves mass. While this may be useful to prevent spurious numerics leading
to crashing care should be taken that the mass doesn't deviate too much.

This scheme is similar to that used by [NEMO-PISCES](https://www.nemo-ocean.eu/), although they scale the
tendency rather than the value, while other Earth system models simply set negative tracers to zero, for
example [NCAR's MARBL](https://marbl-ecosys.github.io/versions/latest_release/index.html) and
[NEMO-TOPAZ2](https://zenodo.org/records/2648099), which does not conserve mass. More complicated schemes
exist, for example [ROMS-BECS](https://zenodo.org/records/3988618) uses an implicite-itterative
approach where each component is updated in sequence to garantee mass conservation, possibly at the
expense of numerical precision.
"""
function ScaleNegativeTracers(tracers; scalefactors = ones(length(tracers)), invalid_fill_value = NaN, warn = false)
    if length(scalefactors) != length(tracers)
        error("Incorrect number of scale factors provided")
    end

    group = NamedTuple{Tuple(tracers)}(Tuple(scalefactors))

    return ScaleNegativeTracers(group; invalid_fill_value, warn)
end

ScaleNegativeTracers(group::NamedTuple{<:Any, <:Tuple{Vararg{Number}}}; invalid_fill_value = NaN, warn = false) =
    ScaleNegativeTracers((; scalefactors = group); invalid_fill_value, warn)

"""
    ScaleNegativeTracers(bgc::AbstractBiogeochemistry; invalid_fill_value = NaN)

Construct a modifier to scale the conserved tracers in `bgc` biogeochemistry.

Groups whose scale factors have both signs are left out: their total (e.g. the oxygen in
`NutrientsPlanktonDetritus` models, which counts organic matter as the oxygen needed to respire it)
can be negative when no tracer is.
"""
function ScaleNegativeTracers(bgc::AbstractBiogeochemistry; invalid_fill_value = NaN, warn = false)
    groups = map(group_scalefactors, conserved_tracers(bgc))

    names = filter(name -> all(≥(0), values(groups[name])), keys(groups))

    return ScaleNegativeTracers(groups[names]; invalid_fill_value, warn)
end

function ScaleNegativeTracers(groups::NamedTuple; invalid_fill_value = NaN, warn = false)
    warn && error("Warning not currently implemented")

    groups = map(group_scalefactors, groups)

    for (group_name, group) in pairs(groups), (name, scalefactor) in pairs(group)
        isfinite(scalefactor) && scalefactor ≥ 0 ||
            throw(ArgumentError("The scale factor of $name in the $group_name group is $scalefactor, but scale factors " *
                                "must not be negative: a total with negative scale factors can be negative when no tracer is."))
    end

    groups = map(group -> (; (name => s for (name, s) in pairs(group) if s > 0)...), groups)
    groups = (; (group_name => group for (group_name, group) in pairs(groups) if !isempty(group))...)

    tracers = Tuple(unique(Symbol[name for group in groups for name in keys(group)]))

    scalefactors = map(group -> map(name -> float(get(group, name, 0)), tracers), groups)

    shared = map(n -> count(s -> s[n] > 0, values(scalefactors)) > 1, Tuple(1:length(tracers)))

    first_own_tracer = map(keys(scalefactors), values(scalefactors)) do group_name, s
        n = findfirst(n -> (s[n] > 0) & !shared[n], 1:length(tracers))

        isnothing(n) && throw(ArgumentError("All of the tracers in the $group_name group are also in other groups, but each " *
                                            "group needs a tracer of its own to make up the changes to its total."))

        return n
    end

    return ScaleNegativeTracers(tracers, scalefactors, shared, first_own_tracer, invalid_fill_value)
end

group_scalefactors(group::NamedTuple) = group
group_scalefactors(names::Tuple{Vararg{Symbol}}) = NamedTuple{names}(map(_ -> 1, names))

summary(scaler::ScaleNegativeTracers) = string("Mass conserving negative scaling of $(scaler.tracers)")

function show(io::IO, scaler::ScaleNegativeTracers)
    print(io, summary(scaler))

    for (n, (group_name, scalefactors)) in enumerate(pairs(scaler.scalefactors))
        group = (; (name => s for (name, s) in zip(scaler.tracers, scalefactors) if s > 0)...)

        print(io, "\n", n == length(scaler.scalefactors) ? "└── " : "├── ", group_name, ": ", group)
    end
end

function update_biogeochemical_state!(model, scale::ScaleNegativeTracers)
    isempty(scale.tracers) && return nothing

    grid = model.grid
    FT = eltype(grid)

    tracers = Tuple(model.tracers[name] for name in scale.tracers)
    scalefactors = map(s -> map(sₙ -> convert(FT, sₙ), s), values(scale.scalefactors))

    launch!(architecture(grid), grid, :xyz, _scale_negative_tracers!,
            map(tracer -> tracer.data, tracers), grid, scalefactors,
            scale.shared, scale.first_own_tracer, convert(FT, scale.invalid_fill_value))

    # `update_state!` filled the halos before calling this, so they still hold the values from before the scaling
    synchronize_communication!(tracers)
    fill_halo_regions!(tracers, model.clock, fields(model))

    return nothing
end

@kernel function _scale_negative_tracers!(tracers, grid, scalefactors, shared, first_own_tracer, invalid_fill_value)
    i, j, k = @index(Global, NTuple)

    scale_negative_tracers!(i, j, k, grid, tracers, scalefactors, shared, first_own_tracer, invalid_fill_value)
end

# The tracers are passed as a tuple and every loop over them is unrolled with `ntuple`,
# which lets the GPU keep them in registers rather than making a local copy of the tuple
@inline function scale_negative_tracers!(i, j, k, grid, tracers::NTuple{M, Any}, args...) where M
    c = ntuple(n -> @inbounds(tracers[n][i, j, k]), Val(M))

    if !immersed_cell(i, j, k, grid) && reduce(|, ntuple(n -> !(c[n] ≥ 0), Val(M)))
        scaled = scale_negative_values(c, args...)

        ntuple(n -> @inbounds(tracers[n][i, j, k] = scaled[n]), Val(M))
    end

    return nothing
end

@inline function scale_negative_values(c::NTuple{M, Any}, s::NTuple{G, Any}, shared,
                                       first_own_tracer, invalid_fill_value) where {M, G}

    c⁺ = ntuple(n -> ifelse(c[n] > 0, c[n], zero(c[n])), Val(M))

    # the total of each group (t), and how much of it is in the positive tracers (p)
    t = ntuple(g -> weighted_sum(s[g], c), Val(G))
    p = ntuple(g -> weighted_sum(s[g], c⁺), Val(G))

    # the factor that would scale the positive tracers of each group to its total if it was the only group
    f = ntuple(g -> ifelse(p[g] > 0, ifelse(t[g] < 0, invalid_fill_value, t[g]) / p[g], one(p[g])), Val(G))

    # a tracer in several groups is scaled by the smallest factor of its groups
    φ = ntuple(n -> reduce(min, ntuple(g -> ifelse(s[g][n] > 0, f[g], one(f[g])), Val(G))), Val(M))

    # the tracers which are only in one group then hold the rest of its total (r): what its positive ones held (q)
    # plus the change (Δ) made by zeroing the negative tracers and scaling the shared ones
    Δ = ntuple(n -> c[n] - ifelse(shared[n], φ[n] * c⁺[n], c⁺[n]), Val(M))
    own = ntuple(n -> ifelse(shared[n], zero(c⁺[n]), c⁺[n]), Val(M))
    scaled_shared = ntuple(n -> ifelse(shared[n], φ[n] * c⁺[n], zero(c⁺[n])), Val(M))

    q = ntuple(g -> weighted_sum(s[g], own), Val(G))
    r = ntuple(g -> ifelse(t[g] < 0, invalid_fill_value - weighted_sum(s[g], scaled_shared),
                                     q[g] + weighted_sum(s[g], Δ)), Val(G))

    # so they are scaled by r / q, or if none of them is positive the first of them gets all of r
    β = ntuple(g -> ifelse(q[g] > 0, r[g] / q[g], zero(r[g])), Val(G))
    remainder = ntuple(g -> ifelse(!(q[g] > 0) & (r[g] > 0), r[g], zero(r[g])) / selected(s[g], first_own_tracer[g]), Val(G))

    return ntuple(Val(M)) do n
        βₙ = sum(ntuple(g -> ifelse(s[g][n] > 0, β[g], zero(β[g])), Val(G)))
        remainderₙ = sum(ntuple(g -> ifelse(first_own_tracer[g] == n, remainder[g], zero(remainder[g])), Val(G)))

        scaled = ifelse(shared[n], φ[n], βₙ) * c[n]
        filled = ifelse(shared[n], zero(c[n]), remainderₙ)

        ifelse(c[n] > 0, scaled, ifelse(isnan(c[n]), c[n], filled))
    end
end

# only the tracers in the group are added, so that e.g. a `NaN` in another group doesn't make the total `NaN`
@inline weighted_sum(s::NTuple{M, Any}, c::NTuple{M, Any}) where M =
    sum(ntuple(n -> ifelse(s[n] > 0, s[n] * c[n], zero(s[n] * c[n])), Val(M)))

# `s[m]` without indexing the tuple with a number that is only known at run time
@inline selected(s::NTuple{M, Any}, m) where M = sum(ntuple(n -> ifelse(n == m, s[n], zero(s[n])), Val(M)))

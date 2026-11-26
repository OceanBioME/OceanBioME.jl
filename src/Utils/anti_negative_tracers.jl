"""
    AntiNegativeTracers(groups; annihilation_tracers, anti_tracer_names)
    AntiNegativeTracers(bgc::AbstractBiogeochemistry; annihilation_tracers, anti_tracer_names)

Correct negative tracer concentrations while conserving elemental inventories by
adding one anti-tracer per conservation group. `groups` maps group names to tracer
stoichiometric weights, following [`ScaleNegativeTracers`](@ref).
`annihilation_tracers` selects an existing tracer in each group whose concentration
can be reduced through annihilation with the anti-tracer.

For example:
```julia
modifier = AntiNegativeTracers(
    (nitrogen = (N = 1, P = 1, Z = 1, D = 1),);
    annihilation_tracers = (nitrogen = :N,))
```

Each anti-tracer is named after its conservation group, followed by `_anti_tracer`
(e.g. `:nitrogen_anti_tracer` or `:phosphate_anti_tracer`). Custom names can be
specified using `anti_tracer_names`. Anti-tracers are automatically added to the
model's tracers when the modifier is used in `Biogeochemistry`.

The correction sets negative tracer concentrations to zero. Any change that
cannot be annihilated using the selected tracer is reflected in the negative
value of the anti-tracer. Annihilation reduces both the concentration of that
tracer and the magnitude of the anti-tracer's negative value. The weighted
inventory is conserved within each grid cell during the correction. Each
anti-tracer uses the same units as the existing tracer selected for annihilation.
Its value is weighted by that tracer's stoichiometric weight when computing the
elemental inventory.

The tracer selected for annihilation must belong exclusively to its conservation
group. Stoichiometric weights are the weights used to calculate elemental
inventories. Groups containing negative weights (such as some oxygen groups)
are skipped by the model-based constructor and rejected by the explicit-group
constructor.

Anti-tracers are ordinary Oceananigans tracers. Configure their transport with
the model's `tracer_advection` and tracer-specific diffusivities in `closure`;
see [Positivity Preservation](@ref pos-preservation) for a runnable example.
"""
struct AntiNegativeTracers{T, R, S, A}
       scalefactors :: S
  annihilation_indices :: A
end

AntiNegativeTracers(tracers::Tuple{Vararg{Symbol}}, scalefactors,
                    anti_tracer_names::Tuple{Vararg{Symbol}}, annihilation_indices) =
    AntiNegativeTracers{tracers, anti_tracer_names, typeof(scalefactors), typeof(annihilation_indices)}(
        scalefactors, annihilation_indices)

@inline Base.getproperty(modifier::AntiNegativeTracers{T, R}, name::Symbol) where {T, R} =
    name === :tracers ? T : name === :anti_tracer_names ? R : getfield(modifier, name)
Base.propertynames(::AntiNegativeTracers, private::Bool=false) =
    (:tracers, :scalefactors, :anti_tracer_names, :annihilation_indices)

function AntiNegativeTracers(bgc::AbstractBiogeochemistry; annihilation_tracers, anti_tracer_names = nothing)
    groups = map(group_scalefactors, conserved_tracers(bgc))
    eligible = filter(name -> all(>=(0), values(groups[name])), keys(groups))

    modifier = AntiNegativeTracers(groups[eligible]; annihilation_tracers, anti_tracer_names)
    existing_tracers = required_biogeochemical_tracers(bgc)
    any(name -> name in existing_tracers, modifier.anti_tracer_names) &&
        throw(ArgumentError("An anti-tracer name conflicts with an existing biogeochemical tracer"))

    return modifier
end

function AntiNegativeTracers(groups::NamedTuple; annihilation_tracers, anti_tracer_names = nothing)
    groups = map(group_scalefactors, groups)
    group_names = keys(groups)

    keys(annihilation_tracers) == group_names ||
        throw(ArgumentError("annihilation_tracers must have exactly the keys $(group_names)"))

    if isnothing(anti_tracer_names)
        anti_tracer_names = NamedTuple{group_names}(map(name -> Symbol(name, :_anti_tracer), group_names))
    end
    keys(anti_tracer_names) == group_names ||
        throw(ArgumentError("anti_tracer_names must have exactly the keys $(group_names)"))

    for (group_name, group) in pairs(groups)
        isempty(group) && throw(ArgumentError("The $group_name group has no tracers"))
        for (name, weight) in pairs(group)
            (isfinite(weight) && weight >= 0) ||
                throw(ArgumentError("The stoichiometric weight of $name in $group_name must be finite and nonnegative"))
        end
        any(>(0), values(group)) ||
            throw(ArgumentError("The $group_name group has no positive stoichiometric weights"))
    end

    tracers = Tuple(unique(Symbol[name for group in values(groups) for (name, s) in pairs(group) if s > 0]))
    weights = map(group -> map(name -> float(get(group, name, 0)), tracers), values(groups))
    names = Tuple(values(anti_tracer_names))

    (all(name -> name isa Symbol, names) && length(unique(names)) == length(names) &&
     all(name -> !(name in tracers), names)) ||
        throw(ArgumentError("Anti-tracer names must be distinct symbols, separate from existing tracer names"))

    annihilation_indices = map(group_names) do group_name
        tracer_name = annihilation_tracers[group_name]
        tracer_name isa Symbol || throw(ArgumentError("The tracer selected for annihilation in $group_name must be a Symbol"))
        n = findfirst(==(tracer_name), tracers)
        isnothing(n) && throw(ArgumentError("The tracer selected for annihilation, $tracer_name, is not in the $group_name group"))
        weights[findfirst(==(group_name), group_names)][n] > 0 ||
            throw(ArgumentError("The tracer selected for annihilation, $tracer_name, has zero weight in $group_name"))

        # Annihilation between the selected tracer and its anti-tracer must not alter
        # any other group (shared tracers themselves are supported).
        count(s -> s[n] > 0, weights) == 1 ||
            throw(ArgumentError("The tracer selected for annihilation, $tracer_name, also belongs to another conservation group"))
        return n
    end

    return AntiNegativeTracers(tracers, weights, names, annihilation_indices)
end

# This method is used by the complete Biogeochemistry wrapper to register
# anti-tracers alongside the underlying biogeochemical tracers.
modifier_tracers(::Any) = ()
modifier_tracers(modifier::AntiNegativeTracers) = modifier.anti_tracer_names
modifier_tracers(modifiers::Tuple) = foldl((names, modifier) -> (names..., modifier_tracers(modifier)...), modifiers; init = ())

summary(modifier::AntiNegativeTracers) =
    "Conservative negative-tracer correction with anti-tracers $(modifier.anti_tracer_names)"

function update_biogeochemical_state!(model, modifier::AntiNegativeTracers)
    isempty(modifier.tracers) && return nothing

    grid = model.grid
    FT = eltype(grid)
    tracers = map(name -> model.tracers[name], modifier.tracers)
    anti_tracers = map(name -> model.tracers[name], modifier.anti_tracer_names)
    weights = map(group -> map(s -> convert(FT, s), group), modifier.scalefactors)

    launch!(architecture(grid), grid, :xyz, _anti_negative_tracers!,
            map(tracer -> tracer.data, tracers), map(anti_tracer -> anti_tracer.data, anti_tracers),
            grid, weights, modifier.annihilation_indices)

    # The model's halos were filled before the modifier was called.
    modified_fields = (tracers..., anti_tracers...)
    synchronize_communication!(modified_fields)
    fill_halo_regions!(modified_fields, model.clock, fields(model))

    return nothing
end

@kernel function _anti_negative_tracers!(tracers, anti_tracers, grid, weights, annihilation_indices)
    i, j, k = @index(Global, NTuple)

    anti_negative_tracers!(i, j, k, grid, tracers, anti_tracers, weights, annihilation_indices)
end

@inline function anti_negative_tracers!(i, j, k, grid, tracers::NTuple{M, Any},
                                          anti_tracers::NTuple{G, Any}, weights, annihilation_indices) where {M, G}
    if !immersed_cell(i, j, k, grid)
        c = ntuple(n -> @inbounds(tracers[n][i, j, k]), Val(M))
        a = ntuple(g -> @inbounds(anti_tracers[g][i, j, k]), Val(G))

        needs_correction = reduce(|, ntuple(n -> !(c[n] >= 0), Val(M))) ||
                           reduce(|, ntuple(g -> a[g] != 0, Val(G)))
        if needs_correction
            corrected, new_anti_tracers = anti_negative_values(c, a, weights, annihilation_indices)
            ntuple(n -> @inbounds(tracers[n][i, j, k] = corrected[n]), Val(M))
            ntuple(g -> @inbounds(anti_tracers[g][i, j, k] = new_anti_tracers[g]), Val(G))
        end
    end
    return nothing
end

# The tuple calculations are statically unrolled so the correction can run in a
# GPU kernel without allocating. Each selected tracer belongs to one group.
@inline function anti_negative_values(c::NTuple{M, Any}, a::NTuple{G, Any},
                                      weights::NTuple{G, Any}, annihilation_indices) where {M, G}
    nonnegative = ntuple(n -> max(c[n], zero(c[n])), Val(M))
    negative_parts = ntuple(n -> c[n] - nonnegative[n], Val(M))
    anti_values = ntuple(g -> a[g] + weighted_sum(weights[g], negative_parts) /
                                  selected(weights[g], annihilation_indices[g]), Val(G))

    # Annihilate each anti-tracer value with its selected nonnegative tracer.
    # Any remaining negative anti-tracer value is conserved.
    combined = ntuple(g -> selected(nonnegative, annihilation_indices[g]) + anti_values[g], Val(G))
    tracer_values = ntuple(g -> max(combined[g], zero(combined[g])), Val(G))
    new_anti_tracers = ntuple(g -> min(combined[g], zero(combined[g])), Val(G))

    corrected = ntuple(Val(M)) do n
        nonnegative[n] + sum(ntuple(g -> ifelse(annihilation_indices[g] == n,
                                                  tracer_values[g] - nonnegative[n],
                                                  zero(nonnegative[n])), Val(G)))
    end

    return corrected, new_anti_tracers
end

"""
    ReservoirNegativeTracers(groups; repayment_tracers, reservoir_names)
    ReservoirNegativeTracers(bgc::AbstractBiogeochemistry; repayment_tracers, reservoir_names)

Construct a conservative negative-tracer correction with one signed reservoir per
conservation group. `groups` maps group names to tracer coefficients, following
[`ScaleNegativeTracers`](@ref). `repayment_tracers` maps each group to a tracer
which can repay its reservoir (typically an inorganic nutrient).

For example:
```julia
modifier = ReservoirNegativeTracers(
    (nitrogen = (N = 1, P = 1, Z = 1, D = 1),);
    repayment_tracers = (nitrogen = :N,))
```

The corresponding reservoir is named `:nitrogen_reservoir` by default. It is
registered as a prognostic tracer when the modifier is used in `Biogeochemistry`.
The biological tracer concentrations are made nonnegative while the weighted
sum of each group's tracers and its reservoir is conserved *in each cell*.

Repayment tracers must be unique to their group: a repayment tracer cannot
contribute to another selected conservation group. Groups with negative
coefficients are not supported; the model-based constructor skips them, as in
`ScaleNegativeTracers(bgc)`.

The elemental-buffer concept follows the negative-buffer approach in the
`ocean_bgc_models` research code.

The reservoirs are ordinary Oceananigans tracers. Configure their transport
with the model's `tracer_advection` and tracer-specific diffusivities in `closure`;
see [Positivity Preservation](@ref pos-preservation) for a runnable example.
"""
struct ReservoirNegativeTracers{T, R, S, A}
       scalefactors :: S
  repayment_indices :: A
end

ReservoirNegativeTracers(tracers::Tuple{Vararg{Symbol}}, scalefactors,
                         reservoir_names::Tuple{Vararg{Symbol}}, repayment_indices) =
    ReservoirNegativeTracers{tracers, reservoir_names, typeof(scalefactors), typeof(repayment_indices)}(
        scalefactors, repayment_indices)

@inline Base.getproperty(modifier::ReservoirNegativeTracers{T, R}, name::Symbol) where {T, R} =
    name === :tracers ? T : name === :reservoir_names ? R : getfield(modifier, name)

function ReservoirNegativeTracers(bgc::AbstractBiogeochemistry; repayment_tracers, reservoir_names = nothing)
    groups = map(group_scalefactors, conserved_tracers(bgc))
    eligible = filter(name -> all(>=(0), values(groups[name])), keys(groups))

    modifier = ReservoirNegativeTracers(groups[eligible]; repayment_tracers, reservoir_names)
    existing_tracers = required_biogeochemical_tracers(bgc)
    any(name -> name in existing_tracers, modifier.reservoir_names) &&
        throw(ArgumentError("A reservoir name conflicts with an existing biogeochemical tracer"))

    return modifier
end

function ReservoirNegativeTracers(groups::NamedTuple; repayment_tracers, reservoir_names = nothing)
    groups = map(group_scalefactors, groups)
    group_names = keys(groups)

    keys(repayment_tracers) == group_names ||
        throw(ArgumentError("repayment_tracers must have exactly the keys $(group_names)"))

    if isnothing(reservoir_names)
        reservoir_names = NamedTuple{group_names}(map(name -> Symbol(name, :_reservoir), group_names))
    end
    keys(reservoir_names) == group_names ||
        throw(ArgumentError("reservoir_names must have exactly the keys $(group_names)"))

    for (group_name, group) in pairs(groups)
        isempty(group) && throw(ArgumentError("The $group_name group has no tracers"))
        for (name, coefficient) in pairs(group)
            (isfinite(coefficient) && coefficient >= 0) ||
                throw(ArgumentError("The coefficient of $name in $group_name must be finite and nonnegative"))
        end
        any(>(0), values(group)) ||
            throw(ArgumentError("The $group_name group has no positive coefficients"))
    end

    tracers = Tuple(unique(Symbol[name for group in values(groups) for (name, s) in pairs(group) if s > 0]))
    weights = map(group -> map(name -> float(get(group, name, 0)), tracers), values(groups))
    names = Tuple(values(reservoir_names))

    (all(name -> name isa Symbol, names) && length(unique(names)) == length(names) &&
     all(name -> !(name in tracers), names)) ||
        throw(ArgumentError("Reservoir names must be distinct symbols, separate from biological tracer names"))

    repayment_indices = map(group_names) do group_name
        anchor = repayment_tracers[group_name]
        anchor isa Symbol || throw(ArgumentError("The repayment tracer for $group_name must be a Symbol"))
        n = findfirst(==(anchor), tracers)
        isnothing(n) && throw(ArgumentError("The repayment tracer $anchor is not in the $group_name group"))
        weights[findfirst(==(group_name), group_names)][n] > 0 ||
            throw(ArgumentError("The repayment tracer $anchor has zero weight in $group_name"))

        # Transferring mass between an anchor and its reservoir must not alter
        # any other group (shared biological tracers themselves are supported).
        count(s -> s[n] > 0, weights) == 1 ||
            throw(ArgumentError("The repayment tracer $anchor also belongs to another conservation group"))
        return n
    end

    return ReservoirNegativeTracers(tracers, weights, names, repayment_indices)
end

# This method is used by the complete Biogeochemistry wrapper to register
# reservoir fields alongside the underlying biogeochemical tracers.
modifier_tracers(::Any) = ()
modifier_tracers(modifier::ReservoirNegativeTracers) = modifier.reservoir_names
modifier_tracers(modifiers::Tuple) = foldl((names, modifier) -> (names..., modifier_tracers(modifier)...), modifiers; init = ())

summary(modifier::ReservoirNegativeTracers) =
    "Conservative negative-tracer correction with reservoirs $(modifier.reservoir_names)"

function update_biogeochemical_state!(model, modifier::ReservoirNegativeTracers)
    isempty(modifier.tracers) && return nothing

    grid = model.grid
    FT = eltype(grid)
    tracers = map(name -> model.tracers[name], modifier.tracers)
    reservoirs = map(name -> model.tracers[name], modifier.reservoir_names)
    weights = map(group -> map(s -> convert(FT, s), group), modifier.scalefactors)

    launch!(architecture(grid), grid, :xyz, _reservoir_negative_tracers!,
            map(tracer -> tracer.data, tracers), map(reservoir -> reservoir.data, reservoirs),
            grid, weights, modifier.repayment_indices)

    # The model's halos were filled before the modifier was called.
    modified_fields = (tracers..., reservoirs...)
    synchronize_communication!(modified_fields)
    fill_halo_regions!(modified_fields, model.clock, fields(model))

    return nothing
end

@kernel function _reservoir_negative_tracers!(tracers, reservoirs, grid, weights, repayment_indices)
    i, j, k = @index(Global, NTuple)

    reservoir_negative_tracers!(i, j, k, grid, tracers, reservoirs, weights, repayment_indices)
end

@inline function reservoir_negative_tracers!(i, j, k, grid, tracers::NTuple{M, Any},
                                           reservoirs::NTuple{G, Any}, weights, repayment_indices) where {M, G}
    if !immersed_cell(i, j, k, grid)
        c = ntuple(n -> @inbounds(tracers[n][i, j, k]), Val(M))
        r = ntuple(g -> @inbounds(reservoirs[g][i, j, k]), Val(G))

        needs_correction = reduce(|, ntuple(n -> !(c[n] >= 0), Val(M))) ||
                           reduce(|, ntuple(g -> r[g] != 0, Val(G)))
        if needs_correction
            corrected, new_reservoirs = reservoir_negative_values(c, r, weights, repayment_indices)
            ntuple(n -> @inbounds(tracers[n][i, j, k] = corrected[n]), Val(M))
            ntuple(g -> @inbounds(reservoirs[g][i, j, k] = new_reservoirs[g]), Val(G))
        end
    end
    return nothing
end

# The tuple calculations are statically unrolled so the correction can run in a
# GPU kernel without allocating. Each repayment tracer is exclusive to one group.
@inline function reservoir_negative_values(c::NTuple{M, Any}, r::NTuple{G, Any},
                                         weights::NTuple{G, Any}, anchors) where {M, G}
    nonnegative = ntuple(n -> max(c[n], zero(c[n])), Val(M))
    deficit = ntuple(n -> c[n] - nonnegative[n], Val(M))
    remaining = ntuple(g -> r[g] + weighted_sum(weights[g], deficit), Val(G))

    # Fold each elemental reservoir back into its designated nonnegative tracer.
    # Any shortfall remains as a negative reservoir, rather than creating mass.
    combined = ntuple(g -> selected(nonnegative, anchors[g]) +
                            remaining[g] / selected(weights[g], anchors[g]), Val(G))
    anchor_values = ntuple(g -> max(combined[g], zero(combined[g])), Val(G))
    new_reservoirs = ntuple(g -> min(combined[g], zero(combined[g])) *
                                 selected(weights[g], anchors[g]), Val(G))

    corrected = ntuple(Val(M)) do n
        nonnegative[n] + sum(ntuple(g -> ifelse(anchors[g] == n,
                                                  anchor_values[g] - nonnegative[n],
                                                  zero(nonnegative[n])), Val(G)))
    end

    return corrected, new_reservoirs
end

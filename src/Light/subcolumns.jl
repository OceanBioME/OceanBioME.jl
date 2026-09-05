#####
##### PAR sub columns
#####
##### A cell's light may be split between sub columns — for example open water and water shaded by ice —
##### each with an area fraction φⱼ and a shortwave ratio rⱼ. Attenuation is shared between them, so
##### PARⱼ(z) = PAR(z)·rⱼ, and any function of light which is *not* linear must therefore be evaluated
##### per sub column and area weighted, Σⱼ φⱼ·f(PAR·rⱼ), rather than applied to the mean light.
#####
##### A light model which splits its light returns a [`SubcolumnPAR`](@ref) as its `PAR` auxiliary field
##### instead of a plain field. The weights therefore travel with the light they weight, and cannot get
##### out of step with it. Everything which just wants the mean light keeps indexing `PAR` as before.
#####
##### For now this is only used in MARBL testing, but in the future we should figure out how to interface
##### the PAR models with it, probably in the NumericalEarth extension
#####

"""
    SubcolumnPAR(field, areas, fractions)

Photosynthetically active radiation split between sub columns — the usual representation of partial ice
cover — with an area fraction φⱼ and a shortwave ratio rⱼ per sub column.

Indexes exactly like the underlying field, returning the mean light, so any rate which is linear in
light needs no changes. Rates which are *not* linear must instead weight over the sub columns with
[`subcolumn_sum`](@ref).

The number of sub columns is uniform over the domain, but the weights need not be — ice cover varies
horizontally — so `areas` and `fractions` are tuples with one entry per sub column, each either a two
dimensional field or a number.
"""
struct SubcolumnPAR{P, A, F}
        field :: P
        areas :: A
    fractions :: F
end

@inline Base.getindex(par::SubcolumnPAR, i, j, k) = @inbounds par.field[i, j, k]
@inline Base.eltype(::SubcolumnPAR{P}) where P = eltype(P)

Adapt.adapt_structure(to, par::SubcolumnPAR) =
    SubcolumnPAR(adapt(to, par.field), adapt(to, par.areas), adapt(to, par.fractions))

summary(par::SubcolumnPAR) = "SubcolumnPAR ($(length(par.areas)) sub columns)"

struct LazyReference{N, V, I}
    values :: V
   indices :: I

    LazyReference(values::NTuple{N},
                  indices::I) where {N, I} =
        new{N, typeof(values), I}(values, indices)
end

@inline Base.getindex(lr::LazyReference, n) = 
    subcolumn_property(lr.values[n], lr.indices)

@inline subcolumn_property(x, indices) = x[indices...]
@inline subcolumn_property(x::Number, indices) = x

struct SubcolumnValues{N, V, A, F}
        value :: V
        areas :: A
    fractions :: F

    function SubcolumnValues(values::V,
                             areas::LazyReference{N},
                             fractions::LazyReference{N}) where {V, N}
        A = typeof(areas)
        F = typeof(fractions)

        return new{N, V, A, F}(values, areas, fractions)
    end
end

# light which is not split just indexes: `@preserve_subcolumns` then costs nothing and the rate sees a
# plain value, so its unsplit method is the one which runs
@inline subcolumn_index(x, i, j, k) = @inbounds x[i, j, k]

@inline subcolumn_index(par::SubcolumnPAR, i, j, k) = 
    @inbounds SubcolumnValues(par.field[i, j, k],
                              LazyReference(par.areas, (i, j, k)),
                              LazyReference(par.fractions, (i, j, k)))

@inline Base.getindex(sv::SubcolumnValues, n) =
    @inbounds sv.value * sv.fractions[n]

@inline Base.getindex(sv::SubcolumnValues, i, j, k) =
    @inbounds sv.value
                              

macro preserve_subcolumns(ex)
    ex isa Expr && ex.head == :ref ||
        error("@preserve_subcolumns expects an indexing expression")

    x = ex.args[1]
    inds = ex.args[2:end]

    # splice the function itself, so a call site needs no import for the macro to work
    return :($subcolumn_index($(esc(x)), $(esc.(inds)...)))
end

# a rate may depend on more than one light value — the two interfaces bounding a cell, say — in which
# case every one of them is evaluated in the same sub column before being weighted
@inline subcolumn_sum(f::F, PAR::SubcolumnValues{N}, others::SubcolumnValues{N}...) where {F, N} =
    sum(ntuple(n -> PAR.areas[n] * f(PAR.value * PAR.fractions[n],
                                     map(o -> o.value * o.fractions[n], others)...), Val(N)))

macro subcolumn_average(def)
    inner = def
    # find the function definition inside other macros
    while Meta.isexpr(inner, :macrocall) 
        inner = last(inner.args)
    end

    Meta.isexpr(inner, (:function, :(=)), 2) && Meta.isexpr(inner.args[1], :call) ||
        error("@subcolumn_average expects a function definition, got: $inner")

    # every argument named `PAR…` carries sub columns and is evaluated in each of them; the rest are
    # passed through unchanged
    call = inner.args[1]
    f, args = call.args[1], call.args[2:end]

    name_of(a) = a isa Symbol ? a : a.args[1]
    is_par(a) = startswith(String(name_of(a)), "PAR")

    par_positions = findall(is_par, args)

    isempty(par_positions) &&
        error("@subcolumn_average expects at least one argument named `PAR…`, got: $call")

    # splice the type and function objects rather than their names, so a call site needs no
    # import of the sub column machinery for the generated method to resolve
    sig = copy(args)

    for p in par_positions
        sig[p] = :($(name_of(args[p]))::$SubcolumnValues)
    end

    # inside the wrapper each sub column value arrives as its own local
    cal = map(name_of, args)
    locals = [Symbol(name_of(args[p]), :ₙ) for p in par_positions]

    for (p, l) in zip(par_positions, locals)
        cal[p] = l
    end

    weighted = Expr(:call, subcolumn_sum,
                    Expr(:->, Expr(:tuple, locals...), Expr(:call, f, cal...)),
                    map(p -> name_of(args[p]), par_positions)...)

    wrapper = Expr(:function, Expr(:call, f, sig...), quote
                       @inline
                       $weighted
                   end)

    return esc(Expr(:block, def, wrapper))
end

#####
##### Interface light
#####

# some rates need the light at the interface between two cells rather than at the cell centre. A light
# model may supply it directly; otherwise it is interpolated from the cell centres.
@inline interface_par(i, j, k, grid, aux::NamedTuple{names}) where {names} =
    _interface_par(i, j, k, grid, aux, aux.PAR, Val(:PAR_interface in names))

@inline _interface_par(i, j, k, grid, aux, PAR, ::Val{true}) = @inbounds aux.PAR_interface[i, j, k]
@inline _interface_par(i, j, k, grid, aux, PAR, ::Val{false}) = ℑzᵃᵃᶠ(i, j, k, grid, PAR)

# the interface light belongs to the same column as the light it bounds, so it is split the same way. A
# supplied `PAR_interface` field is the MEAN interface light, exactly as `PAR` is the mean, so it takes
# its weights from `PAR` rather than carrying its own.
@inline _interface_par(i, j, k, grid, aux, PAR::SubcolumnPAR, ::Val{true}) =
    subcolumn_values(@inbounds(aux.PAR_interface[i, j, k]), PAR, i, j, k)

@inline _interface_par(i, j, k, grid, aux, PAR::SubcolumnPAR, ::Val{false}) =
    subcolumn_values(ℑzᵃᵃᶠ(i, j, k, grid, PAR), PAR, i, j, k)

@inline subcolumn_values(value, PAR::SubcolumnPAR, i, j, k) =
    SubcolumnValues(value,
                    LazyReference(PAR.areas, (i, j, k)),
                    LazyReference(PAR.fractions, (i, j, k)))

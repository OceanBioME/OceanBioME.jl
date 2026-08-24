# check that nothing is zero
@inline autotroph_present(val_prefix::Val, i, j, k, p::ManyPhytoZoo, fields) =
    _autotroph_present(Val(traits(p, val_prefix).silicifier), val_prefix, i, j, k, fields)

@inline _autotroph_present(::Val{false}, val_prefix, i, j, k, fields) =
    @inbounds (autotroph_carbon(val_prefix, fields)[i, j, k]      != 0) &
              (autotroph_chlorophyll(val_prefix, fields)[i, j, k] != 0) &
              (autotroph_phosphorus(val_prefix, fields)[i, j, k]  != 0) &
              (autotroph_iron(val_prefix, fields)[i, j, k]        != 0)

@inline _autotroph_present(::Val{true}, val_prefix, i, j, k, fields) =
    @inbounds _autotroph_present(Val(false), val_prefix, i, j, k, fields) &
              (autotroph_silicon(val_prefix, fields)[i, j, k] != 0)

# get the concentrations of each field
for thing in (:carbon, :chlorophyll, :phosphorus, :iron, :silicon, :calcite)
    func = Symbol(:autotroph_, thing)
    name_fn = Symbol(thing, :_name)
    @eval @generated $func(::Val{S}, fields) where S = :(fields.$($name_fn(S)))
    @eval @inline $func(val_name, i, j, k, p::ManyPhytoZoo, fields) =
        @inbounds $func(val_name, fields)[i, j, k] * autotroph_present(val_name, i, j, k, p, fields)

 end
 
@generated zooplankton_carbon(::Val{Z}, fields) where Z = :(fields.$(carbon_name(Z)))
@inline zooplankton_carbon(val_name, i, j, k, p, fields) =
    @inbounds zooplankton_carbon(val_name, fields)[i, j, k]

#####
##### Generation-time helpers
#####
##### These operate on argument TYPES, not values, and are used by the `@generated` reductions and the
##### tracer list — so they must be defined before any `@generated` function that calls them (world age).
##### The names and traits are all type parameters, so every reduction below unrolls at compile time and
##### an Nₐ×N_z config costs no dispatch. Traits come off each `Autotroph`'s own type, so the same name
##### may calcify in one config and not in another.
#####
##### `phytoplankton_base_names`, `zooplankton_base_names` and `all_autotroph_traits` live in
##### `manifest_many_phyto_zoo.jl`; everything here builds on them.
#####

# reach through a bgc type to the plankton type, so a reduction can take either
plankton_base_type(T) = T <: ManyPhytoZoo ? T : T.parameters[3]

# the traits of one autotroph, by name rather than by position
function base_name_traits(T, s)
    P = plankton_base_type(T)
    traits = all_autotroph_traits(P)

    for (n, name) in phytoplankton_base_names(P)
        name === s && return traits[n].parameters[1]
    end

    error("$s is not an autotroph of this ManyPhytoZoo")
end

phytoplankton_names(T) = Tuple(s for (_, s) in phytoplankton_base_names(plankton_base_type(T)))
zooplankton_names(T)   = Tuple(z for (_, z) in zooplankton_base_names(plankton_base_type(T)))

calcifier_base_names(T)  = Tuple(s for s in phytoplankton_names(T)
                                 if (t = base_name_traits(T, s); t.implicit_calcifier || t.explicit_calcifier))
silicifier_base_names(T) = Tuple(s for s in phytoplankton_names(T) if base_name_traits(T, s).silicifier)
nitrogen_fixer_base_names(T)     = Tuple(s for s in phytoplankton_names(T) if base_name_traits(T, s).nitrogen_fixer)

# the `Heterotroph` type of one predator, and hence its prey classes
function heterotroph_base_type(T, z)
    P = plankton_base_type(T)

    for (n, name) in zooplankton_base_names(P)
        name === z && return P.parameters[5].parameters[2].parameters[n]
    end

    error("$z is not a zooplankton of this ManyPhytoZoo")
end

# the prey-class keys of predator `z`, and the prey members of one class
prey_class_base_names(T, z) = heterotroph_base_type(T, z).parameters[1].parameters[1]

function prey_class_base_type(T, z, c)
    GR = heterotroph_base_type(T, z).parameters[1]

    return GR.parameters[2].parameters[findfirst(==(c), prey_class_base_names(T, z))]
end

prey_base_names(T, z, c) = prey_class_base_type(T, z, c).parameters[1]

# every (predator, prey class) pair, and the ones whose class contains prey `m`
grazing_relationships(T) =
    [(z, c) for z in zooplankton_names(T) for c in prey_class_base_names(T, z)]

grazing_relationships_on(T, m) =
    [(z, c) for (z, c) in grazing_relationships(T) if m in prey_base_names(T, z, c)]

#####
##### Reductions over the members
#####

@generated function sum_over_phytoplankton(f, i, j, k, grid, bgc, fields, aux)
    exprs = [:(f(Val($(QuoteNode(s))), i, j, k, grid, bgc, fields, aux)) for s in phytoplankton_names(bgc)]
    return Expr(:call, :+, exprs...)
end

# same reduction for the plankton-signature helpers (no bgc/aux)
@generated function sum_over_phytoplankton(f, i, j, k, grid, p::ManyPhytoZoo, fields)
    exprs = [:(f(Val($(QuoteNode(s))), i, j, k, grid, p, fields)) for s in phytoplankton_names(p)]
    return Expr(:call, :+, exprs...)
end

@generated function sum_over_zooplankton(f, i, j, k, grid, p::ManyPhytoZoo, fields)
    exprs = [:(f(Val($(QuoteNode(z))), i, j, k, grid, p, fields)) for z in zooplankton_names(p)]
    return Expr(:call, :+, exprs...)
end

# reductions over a trait-selected subset (empty → 0 of the model float type)
for (reduction, names) in ((:sum_over_calcifiers,  :calcifier_base_names),
                           (:sum_over_silicifiers, :silicifier_base_names),
                           (:sum_over_nitrogen_fixers,     :nitrogen_fixer_base_names))
    @eval begin
        @generated function $reduction(f, i, j, k, grid, bgc, fields, aux)
            names = $names(bgc)
            isempty(names) && return :(zero(bgc.plankton.nitrogen_to_carbon))
            exprs = [:(f(Val($(QuoteNode(s))), i, j, k, grid, bgc, fields, aux)) for s in names]
            return Expr(:call, :+, exprs...)
        end

        @generated function $reduction(f, i, j, k, grid, p::ManyPhytoZoo, fields)
            names = $names(p)
            isempty(names) && return :(zero(p.nitrogen_to_carbon))
            exprs = [:(f(Val($(QuoteNode(s))), i, j, k, grid, p, fields)) for s in names]
            return Expr(:call, :+, exprs...)
        end
    end
end

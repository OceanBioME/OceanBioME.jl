#####
##### The grazing matrix
#####
##### Grazing is a triple loop: for each PREDATOR, for each PREY CLASS, for each MEMBER of that class.
##### A prey class is a group — autotrophs and/or zooplankton — grazed as one pool of biomass
##### B = Σ prime(m); the class's flux is then split back over its members in proportion to biomass.
#####

# prime biomass of a prey member: P′ for an autotroph, Z′ for a zooplankton
@generated function prey_active_biomass(::Val{M}, i, j, k, grid, p::ManyPhytoZoo, fields) where M
    M in phytoplankton_names(p)   && return :(loss_active_biomass(Val(M), i, j, k, grid, p, fields))
    M in zooplankton_names(p) && return :(zoo_loss_active_biomass(Val(M), i, j, k, grid, p, fields))
    error("$M is neither an autotroph nor a zooplankton of this ManyPhytoZoo")
end

@inline grazing_relationship(p::ManyPhytoZoo, ::Val{Z}, ::Val{C}) where {Z, C} =
    getproperty(getproperty(p.zooplankton, Z).grazing, C)

# biomass of prey class C of predator Z
@generated function prey_class_biomass(::Val{Z}, ::Val{C}, i, j, k, grid, p::ManyPhytoZoo, fields) where {Z, C}
    exprs = [:(prey_active_biomass(Val($(QuoteNode(m))), i, j, k, grid, p, fields))
             for m in prey_base_names(p, Z, C)]
    return Expr(:call, :+, exprs...)
end

# the class's grazing flux; the predator's temperature factor uses ITS own Q10
@inline function prey_class_grazing(::Val{Z}, ::Val{C}, i, j, k, grid, p::ManyPhytoZoo, fields) where {Z, C}
    g = grazing_relationship(p, Val(Z), Val(C))
    z = getproperty(p.zooplankton, Z)

    @inbounds begin
        T  = fields.T[i, j, k]
        Cz = zooplankton_carbon(Val(Z), fields)[i, j, k]
    end

    B    = prey_class_biomass(Val(Z), Val(C), i, j, k, grid, p, fields)
    Tf_z = temperature_scaling(z.temperature_sensitivity, z.temperature_reference, T)

    return grazing_rate(g.maximum_grazing_rate, g.grazing_half_saturation, g.sigmoidal, Tf_z, Cz, B)
end

# the share of that flux taken from ONE member, (prime_m / B)·G. B = 0 ⇒ 0, and it must be `ifelse`
# rather than ×Bool because prime/B is NaN in the dead branch.
@inline function member_grazing(::Val{Z}, ::Val{C}, ::Val{M}, i, j, k, grid, p::ManyPhytoZoo, fields) where {Z, C, M}
    B     = prey_class_biomass(Val(Z), Val(C), i, j, k, grid, p, fields)
    G     = prey_class_grazing(Val(Z), Val(C), i, j, k, grid, p, fields)
    prime = prey_active_biomass(Val(M), i, j, k, grid, p, fields)

    return ifelse(B > zero(B), prime / B * G, zero(G))
end

# a reduction over every (predator, prey class) relationship that grazes prey M — a prey may be eaten
# by several predators
@generated function sum_over_predators(f, ::Val{M}, i, j, k, grid, p::ManyPhytoZoo, fields) where M
    relationships = grazing_relationships_on(p, M)
    isempty(relationships) && return :(zero(p.nitrogen_to_carbon))
    exprs = [:(f(Val($(QuoteNode(z))), Val($(QuoteNode(c))), Val(M), i, j, k, grid, p, fields))
             for (z, c) in relationships]
    return Expr(:call, :+, exprs...)
end

# total grazing loss of prey M, whether it is an autotroph or a zooplankton
@inline grazing_loss(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    sum_over_predators(member_grazing, val_name, i, j, k, grid, p, fields)

@inline zooplankton_grazing(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    grazing_loss(val_prefix, i, j, k, grid, p, fields)

# the loss of a zooplankton to ITS predators (zero unless some predator's prey class lists it)
@inline zooplankton_grazed(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    grazing_loss(val_name, i, j, k, grid, p, fields)

#####
##### Routing of the grazed flux, per relationship: → predator, → DOC, → POC, → DIC (residual)
#####

@inline _member_graze_to_zooplankton(::Val{Z}, ::Val{C}, ::Val{M}, i, j, k, grid, p, fields) where {Z, C, M} =
    grazing_relationship(p, Val(Z), Val(C)).fraction_to_zooplankton *
    member_grazing(Val(Z), Val(C), Val(M), i, j, k, grid, p, fields)

@inline _member_graze_to_doc(::Val{Z}, ::Val{C}, ::Val{M}, i, j, k, grid, p, fields) where {Z, C, M} =
    grazing_relationship(p, Val(Z), Val(C)).fraction_to_dissolved *
    member_grazing(Val(Z), Val(C), Val(M), i, j, k, grid, p, fields)

# the ballast route is gated on the PREY's implicit_calcifier flag: an implicit calcifier's grazed carbon is
# ballasted by its CaCO₃ (so the relationship's graze_poc is ignored); everything else — including an
# explicit calcifier and any zooplankton prey — uses the relationship's graze_poc.
@inline _member_graze_to_poc(::Val{Z}, ::Val{C}, ::Val{M}, i, j, k, grid, p, fields) where {Z, C, M} =
    _graze_poc_fraction(Val(implicit_calcifier_prey(p, Val(M))),
                        Val(Z), Val(C), Val(M), i, j, k, grid, p, fields) *
    member_grazing(Val(Z), Val(C), Val(M), i, j, k, grid, p, fields)

@inline _graze_poc_fraction(::Val{false}, ::Val{Z}, ::Val{C}, ::Val{M}, i, j, k, grid, p, fields) where {Z, C, M} =
    grazing_relationship(p, Val(Z), Val(C)).fraction_to_particulate

@inline _graze_poc_fraction(::Val{true}, ::Val{Z}, ::Val{C}, ::Val{M}, i, j, k, grid, p, fields) where {Z, C, M} =
    calcifier_graze_poc_fraction(p.calcite_ballast_minimum,
                                 standing_calcite_quota(Val(M), i, j, k, grid, p, fields),
                                 p.small_phyto_poc_factor,
                                 loss_active_biomass(Val(M), i, j, k, grid, p, fields),
                                 p.grazed_small_phyto_poc_limit)

# is prey M an autotroph with implicit_calcifier set? (zooplankton prey never are)
@generated implicit_calcifier_prey(p::ManyPhytoZoo, ::Val{M}) where M =
    :($(M in phytoplankton_names(p) && base_name_traits(p, M).implicit_calcifier))

@inline graze_to_zooplankton(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    sum_over_predators(_member_graze_to_zooplankton, val_name, i, j, k, grid, p, fields)

# the grazed carbon of prey M routed to ONE predator Z — the per-predator breakdown of
# `graze_to_zooplankton`. Only the prey classes of THAT predator which contain M contribute.
@generated function graze_to_predator(::Val{M}, ::Val{Z}, i, j, k, grid, p::ManyPhytoZoo, fields) where {M, Z}
    classes = [c for c in prey_class_base_names(p, Z) if M in prey_base_names(p, Z, c)]
    isempty(classes) && return :(zero(p.nitrogen_to_carbon))
    exprs = [:(_member_graze_to_zooplankton(Val(Z), Val($(QuoteNode(c))), Val(M), i, j, k, grid, p, fields))
             for c in classes]
    return Expr(:call, :+, exprs...)
end

@inline graze_to_doc(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    sum_over_predators(_member_graze_to_doc, val_name, i, j, k, grid, p, fields)

@inline graze_to_poc(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    sum_over_predators(_member_graze_to_poc, val_name, i, j, k, grid, p, fields)

# DIC by residual
@inline graze_to_dic(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    grazing_loss(val_name, i, j, k, grid, p, fields) -
    (graze_to_zooplankton(val_name, i, j, k, grid, p, fields) +
     graze_to_poc(val_name, i, j, k, grid, p, fields) +
     graze_to_doc(val_name, i, j, k, grid, p, fields))

@inline grazing_to_zooplankton(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    graze_to_zooplankton(val_prefix, i, j, k, grid, p, fields)

@inline grazing_to_poc(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    graze_to_poc(val_prefix, i, j, k, grid, p, fields)

@inline grazing_to_doc(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    graze_to_doc(val_prefix, i, j, k, grid, p, fields)

@inline grazing_to_dic(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    graze_to_dic(val_prefix, i, j, k, grid, p, fields)

# loss→POC fraction: a calcifier ballasts its linear loss with its calcite quota (either formula),
# otherwise the PFT's own particulate loss fraction. Trait-dispatched, so the calcite path is never
# taken by a PFT with no CaCO₃ tracer.
@inline loss_to_poc(val_prefix, i, j, k, grid, p::ManyPhytoZoo, fields) =
    _loss_poc_fraction(Val(calcifier(autotroph(p, val_prefix))), val_prefix, i, j, k, grid, p, fields) *
    autotroph_linear_loss(val_prefix, i, j, k, grid, p, fields)

@inline _loss_poc_fraction(::Val{true}, val_prefix, i, j, k, grid, p, fields) =
    standing_calcite_quota(val_prefix, i, j, k, grid, p, fields)

@inline _loss_poc_fraction(::Val{false}, val_prefix, i, j, k, grid, p, fields) =
    autotroph(p, val_prefix).particulate_loss_fraction

#####
##### Zooplankton losses and routing
#####

# loss-active zoo biomass Z′, with the zoo depth ramp (no cold-threshold variant)
@inline function zoo_loss_active_biomass(::Val{Z}, i, j, k, grid, p::ManyPhytoZoo, fields) where Z
    z = getproperty(p.zooplankton, Z)

    @inbounds C = zooplankton_carbon(Val(Z), fields)[i, j, k]

    depth = -znode(i, j, k, grid, Center(), Center(), Center())

    f_thres = clamp((p.zoo_loss_threshold_zero_depth - depth) /
                    (p.zoo_loss_threshold_zero_depth - p.zoo_loss_threshold_full_depth),
                    zero(depth), one(depth))

    return max(C - f_thres * z.loss_threshold_concentration, zero(C))
end

@inline function zooplankton_loss(::Val{Z}, i, j, k, grid, p::ManyPhytoZoo, fields) where Z
    z = getproperty(p.zooplankton, Z)

    @inbounds T = fields.T[i, j, k]

    Z′   = zoo_loss_active_biomass(Val(Z), i, j, k, grid, p, fields)
    Tf_z = temperature_scaling(z.temperature_sensitivity, z.temperature_reference, T)

    return zoo_loss_rate(z.quadratic_mortality_rate, z.linear_mortality_rate,
                         p.zoo_aggregation_exponent, Z′, Tf_z)
end

# food-supply-weighted detrital fraction, PER PREDATOR: the weighted mean of each relationship's
# detritus fraction over everything this predator eats — every member of every one of its prey classes,
# autotroph and zooplankton alike. ε is added per member so the weight stays finite at zero grazing.
@inline _graze_weight(::Val{Z}, ::Val{C}, ::Val{M}, i, j, k, grid, p, fields) where {Z, C, M} =
    member_grazing(Val(Z), Val(C), Val(M), i, j, k, grid, p, fields) +
    p.concentration_regularisation * p.growth_rate_regularisation

@inline _graze_weighted_detritus(::Val{Z}, ::Val{C}, ::Val{M}, i, j, k, grid, p, fields) where {Z, C, M} =
    grazing_relationship(p, Val(Z), Val(C)).detritus_fraction *
    _graze_weight(Val(Z), Val(C), Val(M), i, j, k, grid, p, fields)

# a reduction over every (prey class, member) this predator grazes
@generated function sum_over_prey(f, ::Val{Z}, i, j, k, grid, p::ManyPhytoZoo, fields) where Z
    exprs = [:(f(Val($(QuoteNode(Z))), Val($(QuoteNode(c))), Val($(QuoteNode(m))), i, j, k, grid, p, fields))
             for c in prey_class_base_names(p, Z) for m in prey_base_names(p, Z, c)]
    isempty(exprs) && return :(zero(p.nitrogen_to_carbon))
    return Expr(:call, :+, exprs...)
end

@inline detritus_routing_fraction(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    sum_over_prey(_graze_weighted_detritus, val_name, i, j, k, grid, p, fields) /
    sum_over_prey(_graze_weight, val_name, i, j, k, grid, p, fields)

@inline zoo_loss_to_poc(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    detritus_routing_fraction(val_name, i, j, k, grid, p, fields) *
    zooplankton_loss(val_name, i, j, k, grid, p, fields)

@inline zoo_loss_to_doc(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    (1 - p.labile_fraction) * (1 - detritus_routing_fraction(val_name, i, j, k, grid, p, fields)) *
    zooplankton_loss(val_name, i, j, k, grid, p, fields)

@inline zoo_loss_to_dic(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    p.labile_fraction * (1 - detritus_routing_fraction(val_name, i, j, k, grid, p, fields)) *
    zooplankton_loss(val_name, i, j, k, grid, p, fields)

# food assimilated by predator Z: the retained share of everything it eats — its autotroph prey and any
# zooplankton prey
@inline zooplankton_grazing_gain(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    sum_over_prey(_member_graze_to_zooplankton, val_name, i, j, k, grid, p, fields)

# the carbon a zooplankton LOSES to the rest of the model — mortality plus the part of the grazing on it
# that its predator does not assimilate. Its P/Fe leave at the shared zooplankton quotas, so this is the
# quantity those quotas multiply.
@inline zooplankton_liberated_carbon(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    zooplankton_loss(val_name, i, j, k, grid, p, fields) +
    zooplankton_grazed(val_name, i, j, k, grid, p, fields) -
    graze_to_zooplankton(val_name, i, j, k, grid, p, fields)

# zoo carbon routed to POC / DOC / DIC — mortality plus predation (the latter identically zero unless
# something eats zooplankton)
@inline zooplankton_to_poc(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    zoo_loss_to_poc(val_name, i, j, k, grid, p, fields) +
    graze_to_poc(val_name, i, j, k, grid, p, fields)

@inline zooplankton_to_doc(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    zoo_loss_to_doc(val_name, i, j, k, grid, p, fields) +
    graze_to_doc(val_name, i, j, k, grid, p, fields)

@inline zooplankton_to_dic(val_name, i, j, k, grid, p::ManyPhytoZoo, fields) =
    zoo_loss_to_dic(val_name, i, j, k, grid, p, fields) +
    graze_to_dic(val_name, i, j, k, grid, p, fields)

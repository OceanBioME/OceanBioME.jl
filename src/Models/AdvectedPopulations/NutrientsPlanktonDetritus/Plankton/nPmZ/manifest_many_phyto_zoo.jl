# Generate required field names

phytoplankton_base_names(T) = T <: ManyPhytoZoo ? enumerate(T.parameters[2]) : ()
zooplankton_base_names(T) = T <: ManyPhytoZoo ? enumerate(T.parameters[3]) : ()
all_autotroph_traits(T) = T <: ManyPhytoZoo ? T.parameters[4].parameters[2].parameters : ()

pft_traits(pz::ManyPhytoZoo, s) = traits(pz.autotrophs[s])

# this generator figures out `required_biogeochemical_tracers` at compile time
# so during runtime this just returns the tuple with no overhead
@generated function required_biogeochemical_tracers(p::ManyPhytoZoo)
    tracers = Symbol[]
    base_name_traits = all_autotroph_traits(p)
    for (n, s) in phytoplankton_base_names(p)
        append!(tracers, autotroph_tracer_names(s, base_name_traits[n].parameters[1]))
    end
    for (n, z) in zooplankton_base_names(p)
        push!(tracers, carbon_name(z))
    end
    #push!(tracers, :T)
    return Expr(:tuple, QuoteNode.(tracers)...)
end

required_biogeochemical_auxiliary_fields(::ManyPhytoZoo) = (:PAR, )

# manifest the PFT functions
function manifest_plankton!(autotrophs, zooplankton)
    key = (Tuple((s, traits(a)) for (s, a) in pairs(autotrophs)), keys(zooplankton))
    key in _manifested_marbl_plankton && return nothing
    push!(_manifested_marbl_plankton, key)

    for (sname, autotroph) in pairs(autotrophs)
        flags = traits(autotroph)
        S = QuoteNode(sname)

        for (name, flux) in ((carbon_name(sname),      :photosynthesis_flux),
                             (chlorophyll_name(sname), :photoacclimation_flux),
                             (phosphorus_name(sname),  :phosphorus_uptake_flux),
                             (iron_name(sname),        :iron_uptake_flux))
            @eval @inline (bgc::ManyPhytoZoo_NPD)(i, j, k, grid, ::Val{$(QuoteNode(name))}, clock, fields, aux) =
                autotroph_tendency($flux, quota_of(Val($S), Val($(QuoteNode(name)))), Val($S),
                                   i, j, k, grid, bgc, fields, aux)
        end

        flags.silicifier && @eval begin
            @inline (bgc::ManyPhytoZoo_NPD)(i, j, k, grid, ::Val{$(QuoteNode(silicon_name(sname)))}, clock, fields, aux) =
                autotroph_tendency(silicon_uptake_flux, standing_silicon_quota, Val($S),
                                   i, j, k, grid, bgc, fields, aux)
        end

        (flags.implicit_calcifier || flags.explicit_calcifier) && @eval begin
            @inline (bgc::ManyPhytoZoo_NPD)(i, j, k, grid, ::Val{$(QuoteNode(calcite_name(sname)))}, clock, fields, aux) =
                autotroph_tendency(calcite_formation, standing_calcite_quota, Val($S),
                                   i, j, k, grid, bgc, fields, aux)
        end

        for name in autotroph_tracer_names(sname, flags)
            @eval biogeochemical_drift_velocity(::ManyPhytoZoo_NPD, ::Val{$(QuoteNode(name))}) =
                (u = ZeroField(), v = ZeroField(), w = ZeroField())
        end
    end

    for zname in keys(zooplankton)
        Z = QuoteNode(zname)
        @eval begin
            @inline (bgc::ManyPhytoZoo_NPD)(i, j, k, grid, ::Val{$(QuoteNode(carbon_name(zname)))}, clock, fields, aux) =
                zooplankton_grazing_gain(Val($Z), i, j, k, grid, bgc.plankton, fields) -
                zooplankton_grazed(Val($Z), i, j, k, grid, bgc.plankton, fields) -
                zooplankton_loss(Val($Z), i, j, k, grid, bgc.plankton, fields)

            biogeochemical_drift_velocity(::ManyPhytoZoo_NPD, ::Val{$(QuoteNode(carbon_name(zname)))}) =
                (u = ZeroField(), v = ZeroField(), w = ZeroField())
        end
    end

    return nothing
end

@inline autotroph_tendency(flux, quota, val_base_name, i, j, k, grid, bgc, fields, aux) =
    flux(val_base_name, i, j, k, grid, bgc, fields, aux) -
    quota(val_base_name, i, j, k, grid, bgc.plankton, fields) *
    autotroph_aggregate_loss(val_base_name, i, j, k, grid, bgc.plankton, fields)

@generated function quota_of(::Val{S}, ::Val{name}) where {S, name}
    name === carbon_name(S)      && return :(unit_quota)
    name === chlorophyll_name(S) && return :(chlorophyll_to_carbon_ratio)
    name === phosphorus_name(S)  && return :(standing_phosphorus_quota)
    name === iron_name(S)        && return :(standing_iron_quota)
    error("no quota conjugate to $name")
end

@inline unit_quota(val_base_name, i, j, k, grid, p::ManyPhytoZoo, fields) = one(p.nitrogen_to_carbon)

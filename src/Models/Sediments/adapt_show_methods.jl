Adapt.adapt_structure(to, s::BurialDenitrification) = s   # isbits

Base.summary(::BurialDenitrification) =
    "BurialDenitrification (organic/silica/calcite burial, sedimentary denitrification, anaerobic remineralisation)"

function Base.show(io::IO, s::BurialDenitrification)
    msg  = summary(s) * "\n"
    msg *= "├── organic burial: Dunne et al. (2007), capped at $(s.maximum_organic_burial_fraction)\n"
    msg *= "├── silica burial: Ragueneau et al. (2000), capped at $(s.maximum_silica_burial_fraction)\n"
    msg *= "├── calcite preserved above Ω = $(s.calcite_burial_saturation_threshold)\n"
    msg *= "└── returns: $(coupled_tracers(s))"

    print(io, msg)

    return nothing
end

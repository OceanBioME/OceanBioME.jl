Adapt.adapt_structure(to, nutrients::Nutrients) = 
    Nutrients(adapt(to, nutrients.nitrogen),
              adapt(to, nutrients.phosphate),
              adapt(to, nutrients.iron),
              adapt(to, nutrients.silicate))

Base.summary(nutrients::Nutrients) = string("Nutrients $(required_biogeochemical_tracers(nutrients))")

function Base.show(io::IO, n::Nutrients)
    msg = "Nutrients\n"

    names, values = present_elements(n)
    
    for (n, name) in enumerate(names[1:end-1])
        msg *= "├── $name : $(summary(values[n]))\n"
    end

    msg *= "└── $(names[end]) : $(summary(values[end]))"

    print(io, msg)

    return nothing
end

function present_elements(n::Nutrients)
    vals = (n.nitrogen,
            n.phosphate,
            n.iron,
            n.silicate)

    names = propertynames(n)

    names = [ifelse(isnothing(vals[n]), nothing, name) for (n, name) in enumerate(names)]
    
    names = tuple(unique((nothing, names...))[2:end]...)
    vals = tuple(unique((nothing, vals...))[2:end]...)

    return names, vals
end
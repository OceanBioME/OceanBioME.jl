import Base: show, summary

function Base.show(io::IO, sediment::BiogeochemicalSediment) 

    msg = summary(sediment)

    msg *= "\n    Prognostic fields: $(keys(sediment.fields))"

    msg *= "\n    Tracked fields: $(keys(sediment.tracked_fields))"

    msg *= "\n    Coupled fields: $(coupled_tracers(sediment))"

    print(io, msg)

    return nothing
end

Base.summary(sediment::BiogeochemicalSediment) = 
    string("`BiogeochemicalSediment` with `", summary(sediment.biogeochemistry), "` biogeochemsitry")

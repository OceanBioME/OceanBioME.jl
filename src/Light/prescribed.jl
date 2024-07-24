using Oceananigans.Fields: compute!

maybe_named_fields(field) = ((:PAR, ), (field, ))
maybe_named_fields(fields::NamedTuple) = (keys(fields), values(fields))

struct PrescribedPhotosyntheticallyActiveRadiation{F, FN}
    fields :: F
    field_names :: FN

    PrescribedPhotosyntheticallyActiveRadiation(fields::F, names::FN) where {F, FN} =
        new{F, FN}(fields, names)

    function PrescribedPhotosyntheticallyActiveRadiation(fields)
        names, values = maybe_named_fields(fields)

        F  = typeof(values)
        FN = typeof(names)
        
        return new{F, FN}(values, names)
    end
end

function update_biogeochemical_state!(model, PAR::PrescribedPhotosyntheticallyActiveRadiation)
    for field in PAR.fields
        compute!(field)
    end

    return nothing
end

summary(::PrescribedPhotosyntheticallyActiveRadiation{FT}) where {FT} = string("Prescribed PAR")
show(io::IO, model::PrescribedPhotosyntheticallyActiveRadiation{FT}) where {FT} = print(io, summary(model))

biogeochemical_auxiliary_fields(par::PrescribedPhotosyntheticallyActiveRadiation) = NamedTuple{par.field_names}(par.fields)

adapt_structure(to, par::PrescribedPhotosyntheticallyActiveRadiation) = PrescribedPhotosyntheticallyActiveRadiation(adapt(to, par.fields), 
                                                                                                                    adapt(to, par.field_names))

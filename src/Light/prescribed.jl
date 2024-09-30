using Oceananigans.Architectures: architecture, GPU
using Oceananigans.Fields: compute!, AbstractField, ConstantField

function maybe_named_fields(field)

    isa(field, AbstractField) || @warn "fields: $field is not an `AbstractField"
    
    return NamedTuple{(:PAR, )}((field, ))
end

maybe_named_fields(fields::NamedTuple) = fields


"""
    PrescribedPhotosyntheticallyActiveRadiation(fields)

`PrescribedPhotosyntheticallyActiveRadiation` returns "prescribed" PAR
fields which are user specified, e.g. they may be `FunctionField`s or 
`ConstantField`s.

`fields` may either be an `AbstractField` or a `NamedTuple` of names and 
fields which will be returned in `biogeochemical_auxiliary_fields`, if only
one field is present the field will be named `PAR`.
"""
struct PrescribedPhotosyntheticallyActiveRadiation{F}
    fields :: F

    function PrescribedPhotosyntheticallyActiveRadiation(fields)
        fields = maybe_named_fields(fields)

        F = typeof(fields)

        return new{F}(fields)
    end
end

function update_biogeochemical_state!(model, PAR::PrescribedPhotosyntheticallyActiveRadiation)
    for field in values(PAR.fields)
        compute!(field)
    end

    return nothing
end

summary(::PrescribedPhotosyntheticallyActiveRadiation) = string("Prescribed PAR")
show(io::IO, model::PrescribedPhotosyntheticallyActiveRadiation{F}) where {F} = print(io, summary(model), "\n",
                                                                                            "  Fields:", "\n",
                                                                                            "    └── $(keys(model.fields))")

biogeochemical_auxiliary_fields(par::PrescribedPhotosyntheticallyActiveRadiation) = par.fields

@inline prescribed_field_names(field_names, fields) = field_names
@inline prescribed_field_names(::Nothing, fields::NTuple{N}) where N = tuple(:PAR, ntuple(n -> par_symbol(n), Val(N-1))...)

adapt_structure(to, par::PrescribedPhotosyntheticallyActiveRadiation) = 
    PrescribedPhotosyntheticallyActiveRadiation(adapt(to, par.fields))

using Oceananigans.Architectures: architecture, GPU
using Oceananigans.Fields: compute!, AbstractField

function maybe_named_fields(field)

    isa(field, AbstractField) || @warn "fields: $field is not an `AbstractField"
    
    return ((:PAR, ), (field, ))
end

maybe_named_fields(fields::NamedTuple) = (keys(fields), values(fields))

"""
    PrescribedPhotosyntheticallyActiveRadiation(fields)

`PrescribedPhotosyntheticallyActiveRadiation` returns "prescribed" PAR
fields which are user specified, e.g. they may be `FunctionField`s or 
`ConstantField`s.

`fields` may either be an `AbstractField` or a `NamedTuple` of names and 
fields which will be returned in `biogeochemical_auxiliary_fields`, if only
one field is present the field will be named `PAR`.
"""

struct PrescribedPhotosyntheticallyActiveRadiation{F, FN}
         fields :: F
    field_names :: FN

    PrescribedPhotosyntheticallyActiveRadiation(fields::F, names::FN) where {F, FN} =
        new{F, FN}(fields, names)

    function PrescribedPhotosyntheticallyActiveRadiation(fields)
        names, values = maybe_named_fields(fields)

        isa(architecture(values[1]), GPU) || 
            @warn "On `GPU` prescribed fields will be renamed to `PAR`, `PAR₁`, etc., as symbols can not be passed to the GPU. 
Please make sure they are in this order.
(We're assuming that you're only using this for testing purposes.)"

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

summary(::PrescribedPhotosyntheticallyActiveRadiation) = string("Prescribed PAR")
show(io::IO, model::PrescribedPhotosyntheticallyActiveRadiation{F}) where {F} = print(io, summary(model), "\n",
                                                                                            "  Fields:", "\n",
                                                                                            "    └── $(model.field_names)")

biogeochemical_auxiliary_fields(par::PrescribedPhotosyntheticallyActiveRadiation) = 
    NamedTuple{prescribed_field_names(par.field_names, par.fields)}(par.fields)

@inline prescribed_field_names(field_names, fields) = field_names
@inline prescribed_field_names(::Nothing, fields::NTuple{N}) where N = tuple(:PAR, ntuple(n -> par_symbol(n), Val(N-1))...)

adapt_structure(to, par::PrescribedPhotosyntheticallyActiveRadiation) = 
    PrescribedPhotosyntheticallyActiveRadiation(adapt(to, par.fields), nothing)

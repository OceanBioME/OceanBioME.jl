using Oceananigans.Fields: Field, indices, location

function AuxiliaryFields(fields, grid)
    return NamedTuple(field_name => Field{location(field)...}(grid; indices=indices(field)) for (field_name, field) in pairs(fields))
end
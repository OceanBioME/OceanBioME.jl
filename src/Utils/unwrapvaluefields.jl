"""
Make it superclass of a struct to unwrap fields such that typeof(field) == Val{V}
and make the dot access 's.value_field' return the value V inside Val{V}.

Ideally it would be a trait to avoid having to put it in a type hierarchy, but 
for that we would need to commit gross act of type piracy overriding 
`Base.getproperty`.
"""
abstract type UnwrapValueFields end

unwrap_value(::Val{V}) where {V} = V
unwrap_value(v) = v

Base.getproperty(p::UnwrapValueFields, s::Symbol) = unwrap_value(Base.getfield(p, s))





# check that nothing is zero
@inline autotroph_present(val_name::Val, i, j, k, p::Autotroph, fields) =
    _autotroph_present(Val(traits(p, val_name).silicifier), val_name, i, j, k, fields)

@inline _autotroph_present(::Val{false}, ::Val{S}, i, j, k, fields) where S =
    @inbounds (autotroph_carbon(Val(S), fields)[i, j, k]      != 0) &
              (autotroph_chlorophyll(Val(S), fields)[i, j, k] != 0) &
              (autotroph_phosphorus(Val(S), fields)[i, j, k]  != 0) &
              (autotroph_iron(Val(S), fields)[i, j, k]        != 0)

@inline _autotroph_present(::Val{true}, val_name::Val{S}, i, j, k, fields) where S =
    @inbounds _autotroph_present(Val(false), vs, i, j, k, fields) &
              (autotroph_silicon(Val(S), fields)[i, j, k] != 0)

# get the concentrations of each field
for thing in (:carbon, :chlorophyll, :phosphorus, :iron, :silicon, :calcite)
    func = Symbol(:autotroph_, thing)
    name_fn = Symbol(thing, :_name)
    @eval @generated $func(::Val{S}, fields) where S = :(fields.$($name_fn(S)))
    @eval @inline $func(val_name::Val{S}, i, j, k, p::Autotroph, fields) where S =
        @inbounds $func(val_name, fields)[i, j, k] * autotroph_present(val_name, i, j, k, p, fields)

 end
 
@generated zooplankton_carbon(::Val{Z}, fields) where Z = :(fields.$(carbon_name(Z)))
@inline zooplankton_carbon(val_name, i, j, k, p, field) =
    @inbounds zooplankton_carbon(val_name, fields)[i, j, k]


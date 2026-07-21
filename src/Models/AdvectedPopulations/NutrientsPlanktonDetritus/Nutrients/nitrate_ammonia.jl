
#####
##### Separate nitrogen and ammonia
#####
# this can model how you might make e.g. a complicated iron cycle
"""
    NitrateAmmonia(; nitrification_rate = 5.8e-7)

A nitrogen component for the `nitrogen` slot of [`Nutrients`](@ref) that resolves nitrate (`NO₃`) and
ammonia (`NH₄`) as separate tracers, with ammonia converted to nitrate by nitrification. Plankton may
then take up nitrate and ammonia independently and preferentially.

Keyword Arguments
=================

- `nitrification_rate`: the rate at which ammonia is nitrified to nitrate (1/s)
"""
@kwdef struct NitrateAmmonia{FT}
    nitrification_rate::FT = 5.8e-7 # 1/s
end

Adapt.adapt_structure(to, na::NitrateAmmonia) = 
    NitrateAmmonia(adapt(to, na.nitrification_rate))

Base.summary(::NitrateAmmonia) = string("NitrateAmmonia (:NO₃, :NH₄)")
function Base.show(io::IO, na::NitrateAmmonia)
    msg = summary(na) * "\n"

    msg *= "└── nitrification rate: $(na.nitrification_rate)/s"
    print(io, msg)

    return nothing
end

required_biogeochemical_tracers(::NitrateAmmonia) = (:NO₃, :NH₄)

const NitrateAmmoniaNPD{FT} = NutrientsPlanktonDetritus{FT, <:Nutrients{<:NitrateAmmonia}}

@inline (bgc::NitrateAmmoniaNPD)(i, j, k, grid, val_name::Val{:NH₄}, clock, fields, auxiliary_fields) = (
    inorganic_nitrogen_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + inorganic_nitrogen_waste(i, j, k, grid, bgc.detritus, bgc, fields, auxiliary_fields)
  - nutrient_uptake(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
  - nitrification(i, j, k, grid, bgc.nutrients.nitrogen, fields, auxiliary_fields)
)

@inline (bgc::NitrateAmmoniaNPD)(i, j, k, grid, val_name::Val{:NO₃}, clock, fields, auxiliary_fields) = (
    nitrification(i, j, k, grid, bgc.nutrients.nitrogen, fields, auxiliary_fields)
  - nutrient_uptake(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
)

@inline nitrification(i, j, k, grid, nutrients::NitrateAmmonia, fields, auxiliary_fields) =
    @inbounds fields.NH₄[i, j, k] * nutrients.nitrification_rate

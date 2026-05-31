
#####
##### Separate nitrogen and ammonia
#####
# this can model how you might make e.g. a complicated iron cycle
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
    inorganic_nitrogen_waste(bgc.plankton, bgc, i, j, k, fields, auxiliary_fields)
  + inorganic_nitrogen_waste(bgc.detritus, bgc, i, j, k, fields, auxiliary_fields)
  - nutrient_uptake(bgc.plankton, bgc, i, j, k, val_name, fields, auxiliary_fields)
  - nitrification(bgc.nutrients.nitrogen, i, j, k, fields, auxiliary_fields)
)

@inline (bgc::NitrateAmmoniaNPD)(i, j, k, grid, val_name::Val{:NO₃}, clock, fields, auxiliary_fields) = (
    nitrification(bgc.nutrients.nitrogen, i, j, k, fields, auxiliary_fields)
  - nutrient_uptake(bgc.plankton, bgc, i, j, k, val_name, fields, auxiliary_fields)
)

@inline nitrification(nutrients::NitrateAmmonia, i, j, k, fields, auxiliary_fields) =
    @inbounds fields.NH₄[i, j, k] * nutrients.nitrification_rate

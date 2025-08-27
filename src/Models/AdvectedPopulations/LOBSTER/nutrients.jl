@kwdef struct NitrateAmmonia{FT}
    nitrification_rate::FT = 5.8e-7 # 1/s
end

required_biogeochemical_tracers(::NitrateAmmonia) = (:NO₃, :NH₄)

#####

@kwdef struct NitrateAmmoniaIron{FT}
    nitrification_rate::FT = 5.8e-7 # 1/s
end

required_biogeochemical_tracers(::NitrateAmmoniaIron) = (:NO₃, :NH₄, :Fe)

@inline (lobster::LOBSTER{<:NitrateAmmoniaIron})(i, j, k, grid, val_name::Val{:Fe}, clock, fields, auxiliary_fields) =
    - nutrient_uptake(lobster, i, j, k, val_name, fields, auxiliary_fields)

##### common
const INCLUDES_NITRATE_AMMONIA = Union{NitrateAmmonia, NitrateAmmoniaIron}
const LOBSTER_WITH_NITRATE_AMMONIA = Union{LOBSTER{<:NitrateAmmonia}, LOBSTER{<:NitrateAmmoniaIron}}

@inline nitrifcation(nutrients::INCLUDES_NITRATE_AMMONIA, i, j, k, fields) = @inbounds nutrients.nitrification_rate * fields.NH₄[i, j, k]

@inline (lobster::LOBSTER_WITH_NITRATE_AMMONIA)(i, j, k, grid, val_name::Val{:NO₃}, clock, fields, auxiliary_fields) = (
    nitrifcation(lobster.nutrients, i, j, k, fields) 
  - nutrient_uptake(lobster, i, j, k, val_name, fields, auxiliary_fields)
) # done!

@inline (lobster::LOBSTER_WITH_NITRATE_AMMONIA)(i, j, k, grid, val_name::Val{:NH₄}, clock, fields, auxiliary_fields) = (
    biology_inorganic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)
  + detritus_inorganic_nitrogen_waste(lobster, i, j, k, fields, auxiliary_fields)
  - nitrifcation(lobster.nutrients, i, j, k, fields) 
  - nutrient_uptake(lobster, i, j, k, val_name, fields, auxiliary_fields)
) # done!
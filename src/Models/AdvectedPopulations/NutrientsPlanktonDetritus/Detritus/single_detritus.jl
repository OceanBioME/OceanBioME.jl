using Oceananigans.Units
using Oceananigans.Fields: ZeroField, ConstantField

# Kuhn 2015 "detritus"
struct Detritus{FT, SS}
      remineralisation_rate :: FT
             sinking_speeds :: SS
end

Detritus(FT = Float64;
         remineralisation_rate = 0.1213/day,
         sinking_speeds = (u = ZeroField(), v = ZeroField(), w = ConstantField(-2.7489/day))) =
    Detritus(convert(FT, remineralisation_rate),
             sinking_speeds)

const NPSingleD{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:Detritus}

required_biogeochemical_tracers(::Detritus) = (:D, )
required_biogeochemical_auxiliary_fields(::Detritus) = tuple()

function Detritus(grid::AbstractGrid{FT}; 
                  sinking_speed = 2.7489 / day, # m/s
                  open_bottom = true,
                  kwargs...) where FT

    sinking_speeds = setup_velocity_fields((; D = sinking_speed), grid, open_bottom; three_D = true).D

    return Detritus(FT; sinking_speeds, kwargs...)
end

@inline (bgc::NPSingleD)(i, j, k, grid, val_name::Val{:D}, clock, fields, auxiliary_fields) = (
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  - grazing(bgc.plankton, bgc, val_name, fields, auxiliary_fields) 
  - remineralisation(i, j, k, grid, bgc.detritus, fields, auxiliary_fields)
)

@inline remineralisation(i, j, k, grid, detritus::Detritus, fields, auxiliary_fields) = 
    @inbounds detritus.remineralisation_rate * fields.D[i, j, k]

biogeochemical_drift_velocity(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:Detritus}, ::Val{:D}) = 
    bgc.detritus.sinking_speeds

@inline inorganic_waste(i, j, k, grid, detritus::Detritus, bgc, args...) = 
    remineralisation(detritus, args...)

@inline calcite_dissolution(i, j, k, grid, detritus::Detritus, bgc, fields, auxiliary_fields) = (
    remineralisation(i, j, k, grid, detritus, fields, auxiliary_fields) 
  * carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields) 
  * calcite_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
)

# admin

Adapt.adapt_structure(to, detritus::Detritus) =
    Detritus(remineralisation_rate = adapt(to, detritus.remineralisation_rate),
             sinking_speeds = adapt(to, detritus.sinking_speeds))

Base.summary(::Detritus) = "Detritus (:D)"
function Base.show(io::IO, d::Detritus)
    msg = "Detritus\n"
    msg *= "├── Remineralisation rate: $(d.remineralisation_rate)/s\n"
    msg *= "└── Sinking speed: $(d.sinking_speeds.w)"

    print(io, msg)
    return nothing
end
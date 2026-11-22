using Oceananigans.Units
using Oceananigans.Fields: ZeroField, ConstantField

# Kuhn 2015 "detritus"
"""
    Detritus(grid; sinking_speed = 2.7489/day, dissolution_length = nothing, open_bottom = true,
                   remineralisation_rate = 0.1213/day)

A single-class detritus component (after Kuhn et al., 2015) for the `detritus` slot of a
[`NutrientsPlanktonDetritus`](@ref) model. It adds one sinking detritus tracer `D` which accumulates
plankton waste and grazing residue and is remineralised back to the nutrient pool.

Keyword Arguments
=================

- `grid`: (required) the geometry, needed to configure the sinking-speed field
- `sinking_speed`: the downward sinking speed of detritus (m/s), used unless `dissolution_length`
  is given
- `dissolution_length`: the dissolution length scale of detritus (m), a number or a function of
  `(i, j, k, grid, clock, fields)` such as `DepthDependentDissolutionLength`, enabling
  implicit sinking in place of `sinking_speed` (`D` is then not a tracer)
- `open_bottom`: whether detritus can sink out of the bottom of the domain
- `remineralisation_rate`: the rate at which detritus is remineralised to inorganic nutrients (1/s)
"""
struct Detritus{FT, SK} <: AbstractSinkingDetritus{SK}
      remineralisation_rate :: FT
                    sinking :: SK
end

Detritus(FT = Float64;
         remineralisation_rate = 0.1213/day,
         sinking = ExplicitSinking((u = ZeroField(), v = ZeroField(), w = ConstantField(-2.7489/day)))) =
    Detritus{typeof(convert(FT, remineralisation_rate)), typeof(sinking)}(
        convert(FT, remineralisation_rate), sinking)

const NPSingleD{FT} = NutrientsPlanktonDetritus{FT, <:Any, <:Any, <:Detritus}

function Detritus(grid::AbstractGrid{FT};
                  sinking_speed = 2.7489/day,
                  remineralisation_rate = 0.1213/day,
                  dissolution_length = sinking_speed / remineralisation_rate,
                  implicit_sinking = false,
                  store_flux = false,
                  open_bottom = true) where FT

    if implicit_sinking
        sinking = ImplicitSinking(grid, dissolution_length; tracer_names = (:D,), open_bottom, store_flux)
    else
        sv = setup_velocity_fields((; D = sinking_speed), grid, open_bottom; three_D = true).D
        sinking = ExplicitSinking(sv)
    end

    return Detritus{FT, typeof(sinking)}(convert(FT, remineralisation_rate), sinking)
end

required_biogeochemical_tracers(::Detritus) = (:D, )
required_biogeochemical_tracers(::Detritus{<:Any, <:ImplicitSinking}) = ()

required_biogeochemical_auxiliary_fields(::Detritus) = tuple()

@inline (bgc::NPSingleD)(i, j, k, grid, val_name::Val{:D}, clock, fields, auxiliary_fields) = (
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  - grazing(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
  - remineralisation(i, j, k, grid, bgc.detritus, fields, auxiliary_fields)
)

@inline remineralisation(i, j, k, grid, detritus::Detritus, fields, auxiliary_fields) =
    @inbounds detritus.remineralisation_rate * fields.D[i, j, k]

@inline remineralisation(i, j, k, grid, detritus::Detritus{<:Any, <:ImplicitSinking}, fields, auxiliary_fields) =
    @inbounds detritus.sinking.remineralisation.D[i, j, k]

biogeochemical_drift_velocity(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:Detritus{<:Any, <:ExplicitSinking}}, ::Val{:D}) =
    bgc.detritus.sinking.sinking_speeds

@inline inorganic_waste(i, j, k, grid, detritus::Detritus, bgc, args...) =
    remineralisation(i, j, k, grid, detritus, args...)

@inline calcium_carbonate_dissolution(i, j, k, grid, detritus::Detritus, bgc, fields, auxiliary_fields) = (
    remineralisation(i, j, k, grid, detritus, fields, auxiliary_fields)
  * carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
  * calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
)

@inline implicit_sinking_production(i, j, k, grid, ::Detritus, bgc, fields, aux, ::Val{:D}) =
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux) +
    solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux)

# admin

Adapt.adapt_structure(to, detritus::Detritus) =
    Detritus(remineralisation_rate = adapt(to, detritus.remineralisation_rate),
             sinking = adapt(to, detritus.sinking))

Base.summary(::Detritus{<:Any, <:ExplicitSinking}) = "Detritus (:D)"
Base.summary(::Detritus{<:Any, <:ImplicitSinking}) = "Detritus (implicit sinking)"

function Base.show(io::IO, d::Detritus)
    msg = "Detritus\n"
    msg *= "├── Remineralisation rate: $(d.remineralisation_rate)/s\n"
    msg *= "└── Sinking: $(summary(d.sinking))"
    print(io, msg)
    return nothing
end

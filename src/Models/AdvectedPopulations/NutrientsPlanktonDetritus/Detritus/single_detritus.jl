using Oceananigans.Units
using Oceananigans.Fields: ZeroField, ConstantField, CenterField

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
                  sinking_speed = nothing,
                  dissolution_length = nothing,
                  open_bottom = true,
                  remineralisation_rate = 0.1213/day) where FT

    if dissolution_length !== nothing && sinking_speed !== nothing
        throw(ArgumentError("Cannot specify both `sinking_speed` and `dissolution_length`"))
    end

    if dissolution_length !== nothing
        sinking = ImplicitSinking(grid, dissolution_length; tracer_names = (:D,), open_bottom)
    else
        speed = sinking_speed === nothing ? convert(FT, 2.7489/day) : sinking_speed
        sv = setup_velocity_fields((; D = speed), grid, open_bottom; three_D = true).D
        sinking = ExplicitSinking(sv)
    end

    return Detritus{FT, typeof(sinking)}(convert(FT, remineralisation_rate), sinking)
end

# --- tracers ---

required_biogeochemical_tracers(::Detritus) = (:D, )
required_biogeochemical_tracers(::Detritus{<:Any, <:ImplicitSinking}) = ()

required_biogeochemical_auxiliary_fields(::Detritus) = tuple()

# --- tendencies (explicit sinking only — D is a tracer) ---

@inline (bgc::NPSingleD)(i, j, k, grid, val_name::Val{:D}, clock, fields, auxiliary_fields) = (
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  + solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, auxiliary_fields)
  - grazing(i, j, k, grid, val_name, bgc.plankton, bgc, fields, auxiliary_fields)
  - remineralisation(i, j, k, grid, bgc.detritus, fields, auxiliary_fields)
)

# --- remineralisation ---

@inline remineralisation(i, j, k, grid, detritus::Detritus, fields, auxiliary_fields) =
    @inbounds detritus.remineralisation_rate * fields.D[i, j, k]

@inline remineralisation(i, j, k, grid, detritus::Detritus{<:Any, <:ImplicitSinking}, fields, auxiliary_fields) =
    @inbounds detritus.sinking.remineralisation.D[i, j, k]

# --- drift velocity (explicit only) ---

biogeochemical_drift_velocity(bgc::NutrientsPlanktonDetritus{<:Any, <:Any, <:Any, <:Detritus{<:Any, <:ExplicitSinking}}, ::Val{:D}) =
    bgc.detritus.sinking.sinking_speeds

# --- inorganic waste (nutrient remineralisation source) ---

@inline inorganic_waste(i, j, k, grid, detritus::Detritus, bgc, args...) =
    remineralisation(i, j, k, grid, detritus, args...)

@inline calcium_carbonate_dissolution(i, j, k, grid, detritus::Detritus, bgc, fields, auxiliary_fields) = (
    remineralisation(i, j, k, grid, detritus, fields, auxiliary_fields)
  * carbon_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
  * calcium_carbonate_rain_ratio(i, j, k, grid, bgc.plankton, bgc, fields)
)

# --- implicit sinking production ---

@inline implicit_sinking_production(i, j, k, grid, ::Detritus, bgc, fields, aux, ::Val{:D}) =
    dissolved_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux) +
    solid_waste(i, j, k, grid, bgc.plankton, bgc, fields, aux)

# --- admin ---

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

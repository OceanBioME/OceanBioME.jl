using OceanBioME: Biogeochemistry, ScaleNegativeTracers
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

struct NutrientsPlanktonDetritus{FT, NUT, PLA, DET, CAR, OXY} <: AbstractBiogeochemistry
         nutrients :: NUT 
          plankton :: PLA 
          detritus :: DET 
  inorganic_carbon :: CAR 
            oxygen :: OXY

    NutrientsPlanktonDetritus{FT}(nutrients::NUT, plankton::PLA,
                                  detritus::DET, inorganic_carbon::CAR, oxygen::OXY) where {FT, NUT, PLA, DET, CAR, OXY} =
        new{FT, NUT, PLA, DET, CAR, OXY}(nutrients, plankton,
                                         detritus, inorganic_carbon, 
                                         oxygen)
end

const NPD{FT, NUT, PLA, DET, CAR, OXY} = NutrientsPlanktonDetritus{FT, NUT, PLA, DET, CAR, OXY}

# Fallback for tracers without a component-specific method. Components whose tracer names are
# only known at construction time (`DissolvedParticulate` classes, `CarbonateSystem` replicates, …)
# claim theirs via `component_tendency`, resolved from their type parameters at compile time;
# unclaimed tracers (e.g. `T` and `S`) get zero.
@inline (bgc::NPD)(i, j, k, grid, val_name::Val, clock, fields, auxiliary_fields) =
    component_tendency(i, j, k, grid, bgc.nutrients,        val_name, bgc, clock, fields, auxiliary_fields) +
    component_tendency(i, j, k, grid, bgc.plankton,         val_name, bgc, clock, fields, auxiliary_fields) +
    component_tendency(i, j, k, grid, bgc.detritus,         val_name, bgc, clock, fields, auxiliary_fields) +
    component_tendency(i, j, k, grid, bgc.inorganic_carbon, val_name, bgc, clock, fields, auxiliary_fields) +
    component_tendency(i, j, k, grid, bgc.oxygen,           val_name, bgc, clock, fields, auxiliary_fields)

@inline component_tendency(i, j, k, grid, component, val_name, ::NPD{FT}, clock, fields, auxiliary_fields) where FT = zero(FT)

# Likewise for sinking: each component that can sink either owns the tracer and returns its
# velocity, or passes `fallback` on; `nothing` (no sinking) is Oceananigans' default.
@inline biogeochemical_drift_velocity(bgc::NPD, val_name::Val) =
    component_drift_velocity(bgc.detritus, val_name,
                             component_drift_velocity(bgc.inorganic_carbon, val_name, nothing))

@inline component_drift_velocity(component, val_name, fallback) = fallback

required_biogeochemical_tracers(npd::NutrientsPlanktonDetritus) =
    (required_biogeochemical_tracers(npd.nutrients)...,
     required_biogeochemical_tracers(npd.plankton)...,
     required_biogeochemical_tracers(npd.detritus)...,
     required_biogeochemical_tracers(npd.inorganic_carbon)...,
     required_biogeochemical_tracers(npd.oxygen)...)

required_biogeochemical_auxiliary_fields(npd::NutrientsPlanktonDetritus) =
    (required_biogeochemical_auxiliary_fields(npd.nutrients)...,
     required_biogeochemical_auxiliary_fields(npd.plankton)...,
     required_biogeochemical_auxiliary_fields(npd.detritus)...,
     required_biogeochemical_auxiliary_fields(npd.inorganic_carbon)...,
     required_biogeochemical_auxiliary_fields(npd.oxygen)...)

update_biogeochemical_state!(model, component, npd::NutrientsPlanktonDetritus) = nothing

biogeochemical_auxiliary_fields(npd::NutrientsPlanktonDetritus) =
    merge(biogeochemical_auxiliary_fields(npd.nutrients),
          biogeochemical_auxiliary_fields(npd.plankton),
          biogeochemical_auxiliary_fields(npd.detritus),
          biogeochemical_auxiliary_fields(npd.inorganic_carbon),
          biogeochemical_auxiliary_fields(npd.oxygen))

function update_biogeochemical_state!(model, npd::NutrientsPlanktonDetritus)
    update_biogeochemical_state!(model, npd.nutrients, npd)
    update_biogeochemical_state!(model, npd.plankton, npd)
    update_biogeochemical_state!(model, npd.detritus, npd)
    update_biogeochemical_state!(model, npd.inorganic_carbon, npd)
    update_biogeochemical_state!(model, npd.oxygen, npd)

    return nothing
end

Adapt.adapt_structure(to, npd::NutrientsPlanktonDetritus{FT}) where FT =
    NutrientsPlanktonDetritus{FT}(adapt(to, npd.nutrients),
                                  adapt(to, npd.plankton),
                                  adapt(to, npd.detritus),
                                  adapt(to, npd.inorganic_carbon),
                                  adapt(to, npd.oxygen))

Base.summary(npd::NutrientsPlanktonDetritus{FT}) where FT = 
    string("NutrientsPlanktonDetritus{$FT} with $(required_biogeochemical_tracers(npd))")

function show(io::IO, bgc::NutrientsPlanktonDetritus)
    msg = summary(bgc) * "\n"
    msg *= "├── Plankton: $(summary(bgc.plankton))\n"
    msg *= "├── Nutrients: $(summary(bgc.nutrients))\n"

    if isnothing(bgc.inorganic_carbon) & isnothing(bgc.oxygen)
        msg *= "└── "
    else
        msg *= "├── "
    end

    msg *= "Detritus: $(summary(bgc.detritus))\n"

    if isnothing(bgc.inorganic_carbon) & !isnothing(bgc.oxygen)
        msg *= "└── Oxygen: $(summary(bgc.oxygen))"
    elseif isnothing(bgc.oxygen) & !isnothing(bgc.inorganic_carbon)
        msg *= "└── Carbonate system: $(summary(bgc.inorganic_carbon))"
    elseif !isnothing(bgc.inorganic_carbon) & !isnothing(bgc.oxygen)
        msg *= "├── Carbonate system: $(summary(bgc.inorganic_carbon))\n"
        msg *= "└── Oxygen: $(summary(bgc.oxygen))"
    end

    print(io, string(msg))
end

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

@inline (::NPD{FT})(i, j, k, grid, val_name, clock, fields, auxiliary_fields) where FT = zero(FT)

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

import Base: show, summary

##### NutrientsPlanktonDetritus
summary(bgc::NutrientsPlanktonDetritus) = string("$(overall_name(bgc)) model $(required_biogeochemical_tracers(bgc))")

overall_name(::NutrientsPlanktonDetritus) = "NutrientsPlanktonDetritus"
overall_name(::NutrientsPlanktonDetritus{<:NitrateAmmonia, <:PhytoZoo, <:TwoParticleAndDissolved}) = "LOBSTER"
overall_name(::NutrientsPlanktonDetritus{<:Nutrient, <:PhytoZoo, <:Detritus}) = "NPZD"

function show(io::IO, bgc::NutrientsPlanktonDetritus)
    msg = "$(overall_name(bgc)) model\n"
    msg *= "├── Biology: $(summary(bgc.plankton))\n"
    msg *= "├── Nutrients: $(summary(bgc.nutrients))\n"

    if isnothing(bgc.carbonate_system) & isnothing(bgc.oxygen)
        msg *= "└── "
    else
        msg *= "├── "
    end

    msg *= "Detritus: $(summary(bgc.detritus))\n"

    if isnothing(bgc.carbonate_system) & !isnothing(bgc.oxygen)
        msg *= "└── Oxygen: $(summary(bgc.oxygen))"
    elseif isnothing(bgc.oxygen) & !isnothing(bgc.carbonate_system)
        msg *= "└── Carbonate system: $(summary(bgc.carbonate_system))"
    elseif !isnothing(bgc.carbonate_system) & !isnothing(bgc.oxygen)
        msg *= "├── Carbonate system: $(summary(bgc.carbonate_system))\n"
        msg *= "└── Oxygen: $(summary(bgc.oxygen))"
    end

    print(io, string(msg))
end

##### Components
summary(nutrients::NitrateAmmonia) = string("Nitrate and ammonia $(required_biogeochemical_tracers(nutrients))")
summary(nutrients::NitrateAmmoniaIron) = string("Nitrate, ammonia, and iron $(required_biogeochemical_tracers(nutrients))")
summary(nutrients::Nutrient) = string("Nutrient $(required_biogeochemical_tracers(nutrients))")
summary(plankton::PhytoZoo) = string("Phytoplankton and zooplankton $(required_biogeochemical_tracers(plankton))")
summary(detritus::Detritus) = string("Detritus $(required_biogeochemical_tracers(detritus))")
summary(detritus::TwoParticleAndDissolved) = string("Small and large particles, and dissolved organic matter $(required_biogeochemical_tracers(detritus))")
summary(detritus::VariableRedfieldDetritus) = string("Nitrogen and carbon small and large particles, and dissolved organic matter $(required_biogeochemical_tracers(detritus))")
summary(carbonates::CarbonateSystem) = string("Carbonate chemistry $(required_biogeochemical_tracers(carbonates))")
summary(oxygen::Oxygen) = string("Oxygen chemistry $(required_biogeochemical_tracers(oxygen))")

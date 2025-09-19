import Base: show, summary

##### LOBSTER
summary(lobster::LOBSTER) = string("LOBSTER model $(required_biogeochemical_tracers(lobster))")
function show(io::IO, lobster::LOBSTER)
    msg = "LOBSTER model\n"
    msg *= "├── Biology: $(summary(lobster.biology))\n"
    msg *= "├── Nutrients: $(summary(lobster.nutrients))\n"

    if isnothing(lobster.carbonate_system) & isnothing(lobster.oxygen)
        msg *= "└── "
    else
        msg *= "├── "
    end

    msg *= "Detritus: $(summary(lobster.detritus))\n"

    if isnothing(lobster.carbonate_system) & !isnothing(lobster.oxygen)
        msg *= "└── Oxygen: $(summary(lobster.oxygen))"
    elseif isnothing(lobster.oxygen) & !isnothing(lobster.carbonate_system)
        msg *= "└── Carbonate system: $(summary(lobster.carbonate_system))"
    elseif !isnothing(lobster.carbonate_system) & !isnothing(lobster.oxygen)
        msg *= "├── Carbonate system: $(summary(lobster.carbonate_system))\n"
        msg *= "└── Oxygen: $(summary(lobster.oxygen))"
    end

    print(io, string(msg))
end

##### Components
summary(nutrients::NitrateAmmonia) = string("Nitrate and ammonia $(required_biogeochemical_tracers(nutrients))")
summary(nutrients::NitrateAmmoniaIron) = string("Nitrate, ammonia, and iron $(required_biogeochemical_tracers(nutrients))")
summary(biology::PhytoZoo) = string("Phytoplankton and zooplankton $(required_biogeochemical_tracers(biology))")
summary(detritus::TwoParticleAndDissolved) = string("Small and large particles, and dissolved organic matter $(required_biogeochemical_tracers(detritus))")
summary(detritus::VariableRedfieldDetritus) = string("Nitrogen and carbon small and large particles, and dissolved organic matter $(required_biogeochemical_tracers(detritus))")
summary(carbonates::CarbonateSystem) = string("Carbonate chemistry $(required_biogeochemical_tracers(carbonates))")
summary(oxygen::Oxygen) = string("Oxygen chemistry $(required_biogeochemical_tracers(oxygen))")

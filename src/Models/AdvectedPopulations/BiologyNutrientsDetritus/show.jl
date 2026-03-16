import Base: show, summary

##### LOBSTER
summary(bnd::BiologyNutrientsDetritus) = string("LOBSTER model $(required_biogeochemical_tracers(bnd))")
function show(io::IO, bnd::BiologyNutrientsDetritus)
    msg = "LOBSTER model\n"
    msg *= "├── Biology: $(summary(bnd.biology))\n"
    msg *= "├── Nutrients: $(summary(bnd.nutrients))\n"

    if isnothing(bnd.carbonate_system) & isnothing(bnd.oxygen)
        msg *= "└── "
    else
        msg *= "├── "
    end

    msg *= "Detritus: $(summary(bnd.detritus))\n"

    if isnothing(bnd.carbonate_system) & !isnothing(bnd.oxygen)
        msg *= "└── Oxygen: $(summary(bnd.oxygen))"
    elseif isnothing(bnd.oxygen) & !isnothing(bnd.carbonate_system)
        msg *= "└── Carbonate system: $(summary(bnd.carbonate_system))"
    elseif !isnothing(bnd.carbonate_system) & !isnothing(bnd.oxygen)
        msg *= "├── Carbonate system: $(summary(bnd.carbonate_system))\n"
        msg *= "└── Oxygen: $(summary(bnd.oxygen))"
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

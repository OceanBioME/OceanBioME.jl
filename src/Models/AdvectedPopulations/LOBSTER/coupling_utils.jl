@inline function conserved_tracers(::LOBSTER{<:INCLUDES_NITRATE_AMMONIA, <:PhytoZoo, <:TwoParticleAndDissolved}, labeled = false)
    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)

    return labeled ? (; nitrogen) : nitrogen
end

@inline function conserved_tracers(::LOBSTER{<:INCLUDES_NITRATE_AMMONIA, <:PhytoZoo, <:VariableRedfieldDetritus}, labeled = false)
    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)

    return labeled ? (; nitrogen) : nitrogen
end
    
@inline function conserved_tracers(lobster::LOBSTER{<:INCLUDES_NITRATE_AMMONIA, 
                                                    <:PhytoZoo, 
                                                    <:TwoParticleAndDissolved,
                                                    <:CarbonateSystem},
                                   labeled = false)

    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)

    R = lobster.biology.redfield_ratio
    ρ = lobster.biology.carbon_calcate_ratio

    carbon = (tracers = (:P, :Z, :sPOM, :bPOM, :DOM, :DIC),
              scalefactors = ((1 + ρ) * R, R, R, R, R, 1))

    return labeled ? (; nitrogen, carbon) : (nitrogen, carbon)
end
    
@inline function conserved_tracers(lobster::LOBSTER{<:INCLUDES_NITRATE_AMMONIA, 
                                                    <:PhytoZoo, 
                                                    <:VariableRedfieldDetritus,
                                                    <:CarbonateSystem},
                                   labeled = false)

    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)

    R = lobster.biology.redfield_ratio
    ρ = lobster.biology.carbon_calcate_ratio

    carbon = (tracers = (:P, :Z, :sPOC, :bPOC, :DOC, :DIC),
              scalefactors = ((1 + ρ) * R, R, 1, 1, 1, 1))

    return labeled ? (; nitrogen, carbon) : (nitrogen, carbon)
end
    
@inline chlorophyll(bgc::LOBSTER, model) = bgc.biology.phytoplankton_chlorophyll_ratio * model.tracers.P

@inline conserved_tracers(::LOBSTER{<:INCLUDES_NITRATE_AMMONIA, <:PhytoZoo, <:TwoParticleAndDissolved}) = 
    (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)

@inline conserved_tracers(::LOBSTER{<:INCLUDES_NITRATE_AMMONIA, <:PhytoZoo, <:VariableRedfieldDetritus}) = 
    (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)
    
@inline function conserved_tracers(lobster::LOBSTER{<:INCLUDES_NITRATE_AMMONIA, 
                                                    <:PhytoZoo, 
                                                    <:TwoParticleAndDissolved,
                                                    <:CarbonateSystem})

    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)

    R = lobster.biology.redfield_ratio

    carbon = (tracers = (:P, :Z, :sPOM, :bPOM, :DOM, :DIC),
              scalefactors = (R, R, R, R, R, 1))

    return (nitrogen, carbon)
end
    
@inline function conserved_tracers(lobster::LOBSTER{<:INCLUDES_NITRATE_AMMONIA, 
                                                    <:PhytoZoo, 
                                                    <:VariableRedfieldDetritus,
                                                    <:CarbonateSystem})

    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)

    R = lobster.biology.redfield_ratio

    carbon = (tracers = (:P, :Z, :sPOC, :bPOC, :DOC, :DIC),
              scalefactors = (R, R, 1, 1, 1, 1))

    return (nitrogen, carbon)
end
    
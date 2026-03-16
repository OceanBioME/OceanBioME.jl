@inline function conserved_tracers(::BiologyNutrientsDetritus{<:PhytoZoo, 
                                                              <:INCLUDES_NITRATE_AMMONIA, 
                                                              <:TwoParticleAndDissolved}, 
                                   labeled = false)
    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)

    return labeled ? (; nitrogen) : nitrogen
end

@inline function conserved_tracers(::BiologyNutrientsDetritus{<:PhytoZoo, 
                                                              <:INCLUDES_NITRATE_AMMONIA, 
                                                              <:VariableRedfieldDetritus}, 
                                   labeled = false)
    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)

    return labeled ? (; nitrogen) : nitrogen
end
    
@inline function conserved_tracers(bnd::BiologyNutrientsDetritus{<:PhytoZoo,
                                                                 <:INCLUDES_NITRATE_AMMONIA, 
                                                                 <:TwoParticleAndDissolved,
                                                                 <:CarbonateSystem},
                                   labeled = false)

    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPOM, :bPOM, :DOM)

    R = bnd.biology.redfield_ratio
    ρ = bnd.biology.carbon_calcate_ratio

    carbon = (tracers = (:P, :Z, :sPOM, :bPOM, :DOM, :DIC),
              scalefactors = ((1 + ρ) * R, R, R, R, R, 1))

    return labeled ? (; nitrogen, carbon) : (nitrogen, carbon)
end
    
@inline function conserved_tracers(bnd::BiologyNutrientsDetritus{<:PhytoZoo, 
                                                                 <:INCLUDES_NITRATE_AMMONIA,
                                                                 <:VariableRedfieldDetritus,
                                                                 <:CarbonateSystem},
                                   labeled = false)

    nitrogen = (:NO₃, :NH₄, :P, :Z, :sPON, :bPON, :DON)

    R = bnd.biology.redfield_ratio
    ρ = bnd.biology.carbon_calcate_ratio

    carbon = (tracers = (:P, :Z, :sPOC, :bPOC, :DOC, :DIC),
              scalefactors = ((1 + ρ) * R, R, 1, 1, 1, 1))

    return labeled ? (; nitrogen, carbon) : (nitrogen, carbon)
end
    
@inline chlorophyll(bgc::BiologyNutrientsDetritus, model) = bgc.biology.phytoplankton_chlorophyll_ratio * model.tracers.P

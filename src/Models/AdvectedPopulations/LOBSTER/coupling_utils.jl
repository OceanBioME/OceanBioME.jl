@inline function conserved_nitrogen_tracers(lobster::LOBSTER)
    tracers = (:P, :Z)

    if lobster.nutrients isa INCLUDES_NITRATE_AMMONIA
        tracers = (tracers..., :NO₃, :NH₄)
    else
        tracers = (tracers..., :N)
    end

    if lobster.detritus isa VariableRedfieldDetritus
        tracers = (tracers..., :sPON, :bPON, :DON)
    elseif lobster.detritus isa TwoParticleAndDissolved
        tracers = (tracers..., :sPOM, :bPOM, :DOM)
    else
        tracers = (tracers..., :D)
    end

    return tracers
end

@inline conserved_tracers(lobster::LOBSTER, labeled = false) =
    labeled ? (; nitrogen = conserved_nitrogen_tracers(lobster)) : conserved_nitrogen_tracers(lobster)
    
@inline function conserved_tracers(lobster::LOBSTER{<:Any, 
                                                    <:Any, 
                                                    <:Any,
                                                    <:CarbonateSystem},
                                   labeled = false)

    nitrogen = conserved_nitrogen_tracers(lobster)

    R = lobster.biology.redfield_ratio
    ρ = lobster.biology.carbon_calcate_ratio

    carbon_tracers = (:P, :Z, :DIC)
    carbon_scalefactors = ((1 + ρ) * R, R, 1)

    if lobster.detritus isa VariableRedfieldDetritus
        carbon_tracers = (carbon_tracers..., :sPOC, :bPOC, :DOC)
        carbon_scalefactors = (carbon_scalefactors..., 1, 1, 1)
    elseif lobster.detritus isa TwoParticleAndDissolved
        carbon_tracers = (carbon_tracers..., :sPOM, :bPOM, :DOM)
        carbon_scalefactors = (carbon_scalefactors..., R, R, R)
    else
        carbon_tracers = (carbon_tracers..., :D)
        carbon_scalefactors = (carbon_scalefactors..., R)
    end

    carbon = (tracers = carbon_tracers, scalefactors = carbon_scalefactors)

    return labeled ? (; nitrogen, carbon) : (nitrogen, carbon)
end
    
@inline chlorophyll(bgc::LOBSTER, model) = bgc.biology.phytoplankton_chlorophyll_ratio * model.tracers.P

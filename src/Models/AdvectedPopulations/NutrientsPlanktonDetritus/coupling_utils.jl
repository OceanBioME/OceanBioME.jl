@inline function conserved_nitrogen_tracers(bpd::NutrientsPlanktonDetritus)
    tracers = (:P, :Z)

    if bpd.nutrients isa IncludesNitrateAmmonia
        tracers = (tracers..., :NO₃, :NH₄)
    else
        tracers = (tracers..., :N)
    end

    if bpd.detritus isa VariableRedfieldDetritus
        tracers = (tracers..., :sPON, :bPON, :DON)
    elseif bpd.detritus isa TwoParticleAndDissolved
        tracers = (tracers..., :sPOM, :bPOM, :DOM)
    else
        tracers = (tracers..., :D)
    end

    return tracers
end

@inline conserved_tracers(bpd::NutrientsPlanktonDetritus, labeled = false) =
    labeled ? (; nitrogen = conserved_nitrogen_tracers(bpd)) : conserved_nitrogen_tracers(bpd)
    
@inline function conserved_tracers(bpd::NutrientsPlanktonDetritus{<:Any, 
                                                    <:Any, 
                                                    <:Any,
                                                    <:CarbonateSystem},
                                   labeled = false)

    nitrogen = conserved_nitrogen_tracers(bpd)

    R = bpd.plankton.redfield_ratio
    ρ = bpd.plankton.carbon_calcate_ratio

    carbon_tracers = (:P, :Z, :DIC)
    carbon_scalefactors = ((1 + ρ) * R, R, 1)

    if bpd.detritus isa VariableRedfieldDetritus
        carbon_tracers = (carbon_tracers..., :sPOC, :bPOC, :DOC)
        carbon_scalefactors = (carbon_scalefactors..., 1, 1, 1)
    elseif bpd.detritus isa TwoParticleAndDissolved
        carbon_tracers = (carbon_tracers..., :sPOM, :bPOM, :DOM)
        carbon_scalefactors = (carbon_scalefactors..., R, R, R)
    else
        carbon_tracers = (carbon_tracers..., :D)
        carbon_scalefactors = (carbon_scalefactors..., R)
    end

    carbon = (tracers = carbon_tracers, scalefactors = carbon_scalefactors)

    return labeled ? (; nitrogen, carbon) : (nitrogen, carbon)
end
    
@inline chlorophyll(bgc::NutrientsPlanktonDetritus, model) = bgc.plankton.phytoplankton_chlorophyll_ratio * model.tracers.P

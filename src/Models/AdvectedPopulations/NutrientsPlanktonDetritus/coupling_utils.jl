@inline function conserved_nitrogen_tracers(bpd::NutrientsPlanktonDetritus)
    tracers = tuple()

    if bpd.plankton isa PhytoZoo
        tracers = (tracers..., (:P, :Z)...)
    end

    if bpd.nutrients isa IncludesNitrateAmmonia
        tracers = (tracers..., :NO₃, :NH₄)
    else
        tracers = (tracers..., :N)
    end

    if bpd.detritus isa VariableRedfieldDetritus
        tracers = (tracers..., :sPON, :bPON, :DON)
    elseif bpd.detritus isa TwoParticleAndDissolved
        tracers = (tracers..., :sPOM, :bPOM, :DOM)
    elseif bpd.detritus isa Detritus
        tracers = (tracers..., :D)
    end

    return tracers
end

@inline conserved_tracers(bpd::NutrientsPlanktonDetritus, labeled = false) =
    labeled ? (; nitrogen = conserved_nitrogen_tracers(bpd)) : conserved_nitrogen_tracers(bpd)
    
@inline function conserved_tracers(bpd::NutrientsPlanktonDetritus{<:Any, 
                                                    <:Any, 
                                                    <:Any,
                                                    <:CarbonateSystem{N, calcite}},
                                   labeled = false) where {N, calcite}

    nitrogen = conserved_nitrogen_tracers(bpd)

    carbon_tracers = tuple()
    carbon_scalefactors = tuple()

    if bpd.plankton isa PhytoZoo
        R = bpd.plankton.redfield_ratio
        ρ = bpd.plankton.carbon_calcite_ratio

        carbon_tracers = (carbon_tracers..., (:P, :Z)...)
        carbon_scalefactors = (carbon_scalefactors..., (1 + ρ) * R, R)
    else
        R = 6.56
    end

    if bpd.detritus isa VariableRedfieldDetritus
        carbon_tracers = (carbon_tracers..., :sPOC, :bPOC, :DOC)
        carbon_scalefactors = (carbon_scalefactors..., 1, 1, 1)
    elseif bpd.detritus isa TwoParticleAndDissolved
        carbon_tracers = (carbon_tracers..., :sPOM, :bPOM, :DOM)
        carbon_scalefactors = (carbon_scalefactors..., R, R, R)
    elseif bpd.detritus isa Detritus
        carbon_tracers = (carbon_tracers..., :D)
        carbon_scalefactors = (carbon_scalefactors..., R)
    end

    if (N == 1)
        if calcite
            carbon_tracers = (carbon_tracers..., :DIC, :CaCO₃)
            carbon_scalefactors = (carbon_scalefactors..., 1, 1)
        else
            carbon_tracers = (carbon_tracers..., :DIC)
            carbon_scalefactors = (carbon_scalefactors..., 1)
        end
        carbon = (tracers = carbon_tracers, scalefactors = carbon_scalefactors)
    else
        if calcite
            DICs = [Symbol(:DIC, n) for n in 1:N]
            CaCO₃s = [Symbol(:CaCO₃, n) for n in 1:N]
            carbon = [(tracers = (carbon_tracers..., DIC, CaCO₃s[n]), 
                    scalefactors = (carbon_scalefactors..., 1, 1))
                    for (n, DIC) in enumerate(DICs)]
        else
            DICs = [Symbol(:DIC, n) for n in 1:N]
            carbon = [(tracers = (carbon_tracers..., DIC), 
                    scalefactors = (carbon_scalefactors..., 1))
                    for DIC in DICs]
        end
    end

    return labeled ? (; nitrogen, carbon) : (nitrogen, carbon)
end
    
@inline chlorophyll(bgc::NutrientsPlanktonDetritus, model) = bgc.plankton.phytoplankton_chlorophyll_ratio * model.tracers.P

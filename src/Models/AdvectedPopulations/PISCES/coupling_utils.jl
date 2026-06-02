# sediment models
@inline redfield(val_name, bgc::PISCES, tracers) = bgc.nitrogen_redfield_ratio

@inline sinking_tracers(::PISCES) = (:POC, :GOC, :SFe, :BFe, :PSi, :CaCO₃) # please list them here

# light attenuation model
@inline chlorophyll(::PISCES, model) = model.tracers.PChl + model.tracers.DChl

# negative tracer scaling
# TODO: deal with remaining (PChl, DChl, O₂, Alk) - latter two should never be near zero
@inline function conserved_tracers(bgc::PISCES)
    carbon = NamedTuple{(:P, :D, :Z, :M, :DOC, :POC, :GOC, :DIC, :CaCO₃)}(ones(9))

    # iron ratio for DOC might be wrong
    iron = (PFe = 1,
            DFe = 1,
            Z = bgc.zooplankton.micro.iron_ratio,
            M = bgc.zooplankton.meso.iron_ratio,
            SFe = 1,
            BFe = 1,
            Fe = 1)
    

    θ_PO₄ = bgc.phosphate_redfield_ratio
    phosphate = (P = θ_PO₄, D = θ_PO₄, Z = θ_PO₄, M = θ_PO₄,
                 DOC = θ_PO₄, POC = θ_PO₄, GOC = θ_PO₄, PO₄ = 1)

    silicon = NamedTuple{(:DSi, :Si, :PSi)}(ones(3))

    θN = bgc.nitrogen_redfield_ratio
    nitrogen = (NH₄ = 1, NO₃ = 1, P = θN, D = θN, Z = θN, M = θN, DOC = θN, POC = θN, GOC = 1)
    
    return (; carbon, iron, phosphate, silicon, nitrogen)
end
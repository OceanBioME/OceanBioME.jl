# sediment models
@inline redfield(i, j, k, val_tracer_name, bgc::PISCES, tracers) = NaN

@inline nitrogen_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline carbon_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline remineralisation_receiver(::PISCES) = :NH₄

@inline sinking_tracers(::PISCES) = (:POC, :GOC, :SFe, :BFe, :PSi, :CaCO₃) # please list them here

# light attenuation model
@inline chlorophyll(::PISCES, model) = model.tracers.PChl + model.tracers.DChl

# negative tracer scaling
# TODO: make this work in scale negatives (i.e. return multiple functions for each group)
# TODO: deal with remaining (DSi, Si, PChl, DChl, O₂, Alk) - latter two should never be near zero
@inline function conserved_tracers(bgc::PISCES)
    carbon = (:P, :D, :Z, :M, :DOC, :POC, :GOC, :DIC, :CaCO₃)

    # iron ratio for DOC might be wrong
    iron = (tracers = (:PFe, :DFe, :Z, :M, :DOC, :SFe, :BFe, :Fe), 
            scalefactors = (1, 1, bgc.microzooplankton.iron_ratio, bgc.mesozooplankton.iron_ratio, bgc.iron_redfield_ratio, 1, 1, 1))
    
    θ_PO₄ = bgc.phosphate_redfield_ratio
    phosphate = (tracers = (:P, :D, :Z, :M, :DOC, :POC, :GOC, :PO₄),
                 scalefactors = (θ_PO₄, θ_PO₄, θ_PO₄, θ_PO₄, θ_PO₄, θ_PO₄, θ_PO₄, 1))

    silicon = (:DSi, :Si, :PSi)

    return (; carbon, iron, phosphate, silicon)
end
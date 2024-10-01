import OceanBioME.Models.Sediments: nitrogen_flux, carbon_flux, remineralisation_receiver, sinking_tracers

# sediment models
@inline redfield(val_name, bgc::PISCES, tracers) = bgc.nitrogen_redfield_ratio

@inline nitrogen_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = bgc.nitrogen_redfield_ratio * carbon_flux(i, j, k, grid, advection, bgc, tracers)

@inline carbon_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = sinking_flux(i, j, k, grid, adveciton, bgc, Val(:POC), tracers) +
                                                                      sinking_flux(i, j, k, grid, adveciton, bgc, Val(:GOC), tracers)

@inline remineralisation_receiver(::PISCES) = :NH₄

@inline sinking_tracers(::PISCES) = (:POC, :GOC, :SFe, :BFe, :PSi, :CaCO₃) # please list them here

# light attenuation model
@inline chlorophyll(::PISCES, model) = model.tracers.PChl + model.tracers.DChl

# negative tracer scaling
# TODO: deal with remaining (PChl, DChl, O₂, Alk) - latter two should never be near zero
@inline function conserved_tracers(bgc::PISCES; ntuple = false)
    carbon = (:P, :D, :Z, :M, :DOC, :POC, :GOC, :DIC, :CaCO₃)

    # iron ratio for DOC might be wrong
    iron = (tracers = (:PFe, :DFe, :Z, :M, :SFe, :BFe, :Fe), 
            scalefactors = (1, 1, bgc.zooplankton.micro.iron_ratio, bgc.zooplankton.meso.iron_ratio, 1, 1, 1))
    
    θ_PO₄ = bgc.phosphate_redfield_ratio
    phosphate = (tracers = (:P, :D, :Z, :M, :DOC, :POC, :GOC, :PO₄),
                 scalefactors = (θ_PO₄, θ_PO₄, θ_PO₄, θ_PO₄, θ_PO₄, θ_PO₄, θ_PO₄, 1))

    silicon = (:DSi, :Si, :PSi)

    θN = bgc.nitrogen_redfield_ratio
    nitrogen = (tracers = (:NH₄, :NO₃, :P, :D, :Z, :M, :DOC, :POC, :GOC),
                scalefactors = (1, 1, θN, θN, θN, θN, θN, θN, θN))

    if ntuple
        return (; carbon, iron, phosphate, silicon, nitrogen)
    else
        return (carbon, iron, phosphate, silicon, nitrogen)
    end
end
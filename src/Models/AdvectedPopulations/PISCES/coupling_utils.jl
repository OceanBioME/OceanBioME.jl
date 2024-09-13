# sediment models
@inline redfield(i, j, k, val_tracer_name, bgc::PISCES, tracers) = NaN

@inline nitrogen_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline carbon_flux(i, j, k, grid, advection, bgc::PISCES, tracers) = NaN

@inline remineralisation_receiver(::PISCES) = :NH₄

@inline sinking_tracers(::PISCES) = (:POC, :GOC, :SFe, :BFe, :PSi, :CaCO₃) # please list them here

# light attenuation model
@inline chlorophyll(bgc::PISCES, model) = model.tracers.Pᶜʰˡ + model.tracers.Dᶜʰˡ

# negative tracer scaling
# TODO: make this work in scale negatives (i.e. return multiple functions for each group)
# TODO: deal with remaining (DSi, Si, PChl, DChl, O₂, Alk) - latter two should never be near zero
@inline conserved_tracers(::PISCES) = 
    (carbon = (:P, :D, :Z, :M, :DOC, :POC, :GOC, :DIC, :CaCO₃),
     iron = (tracers = (:PFe, :DFe, :Z, :M, DOC, :SFe, :BFe, :Fe), 
             scalefactors = (1, 1, A, A, A, 1, 1, 1)),
     phosphate = (tracers = (:P, :D, :Z, :M, :DOC, :POC, :GOC, :PO₄),
                  scalefactors = (A, A, A, A, A, A, A, 1)))
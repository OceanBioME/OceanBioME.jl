function Adapt.adapt_structure(to, d::MultiElementRefractoryDissolvedParticulate{FT}) where FT
    sinking_velocities = adapt(to, d.sinking_velocities)

    return MultiElementRefractoryDissolvedParticulate{FT,
                                             typeof(d.semilabile_remineralisation_rate_light),
                                             typeof(d.refractory_production_fraction),
                                             typeof(sinking_velocities)}(
        d.semilabile_remineralisation_rate_light,
        d.semilabile_remineralisation_rate_dark,
        d.refractory_remineralisation_rate,
        d.photodegradation_rate,
        d.particulate_organic_remineralisation_rate,
        d.particulate_silica_remineralisation_rate,
        d.refractory_production_fraction,
        d.particulate_refractory_fraction,
        d.uv_reference_depth,
        sinking_velocities)
end

Base.summary(::MultiElementRefractoryDissolvedParticulate) =
    "MultiElementRefractoryDissolvedParticulate (DOC/DON/DOP, DOCr/DONr/DOPr, POC/POP/PFe, bSi)"

function Base.show(io::IO, d::MultiElementRefractoryDissolvedParticulate)
    msg  = summary(d) * "\n"
    msg *= "├── PON is diagnostic (fixed particulate N:C), CaCO₃ lives in the inorganic carbon\n"
    msg *= "├── particulate remineralisation: $(d.particulate_organic_remineralisation_rate) 1/s"
    msg *= " (silica $(d.particulate_silica_remineralisation_rate) 1/s)\n"
    msg *= "└── tracers: $(required_biogeochemical_tracers(d))"

    print(io, msg)

    return nothing
end

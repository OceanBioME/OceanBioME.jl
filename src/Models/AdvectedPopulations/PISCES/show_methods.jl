using OceanBioME.Models.PISCESModel.Phytoplankton: MixedMondo
using OceanBioME.Models.PISCESModel.Zooplankton: QualityDependantZooplankton

summary(bgc::PISCES) = string("PISCES biogeochemical model ($(length(required_biogeochemical_tracers(bgc))-2) tracers)") # hack to exclude temp and salinity

function show(io::IO, bgc::PISCES)

    FT = typeof(bgc.background_shear)   

    output = "Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model {$FT}"

    output *= "\n Phytoplankton: $(summary(bgc.phytoplankton))"

    output *= "\n Zooplankton: $(summary(bgc.zooplankton))"

    output *= "\n Dissolved organic matter: $(summary(bgc.dissolved_organic_matter))"

    output *= "\n Particulate organic matter: $(summary(bgc.particulate_organic_matter))"

    output *= "\n Nitrogen: $(summary(bgc.nitrogen))"

    output *= "\n Iron: $(summary(bgc.iron))"

    output *= "\n Silicate: $(summary(bgc.silicate))"

    output *= "\n Oxygen: $(summary(bgc.oxygen))"

    output *= "\n Phosphate: $(summary(bgc.phosphate))"

    output *= "\n Inorganic carbon: $(summary(bgc.inorganic_carbon))"

    output *= "\n Latitude: $(summary(bgc.latitude))"

    output *= "\n Day length: $(nameof(bgc.day_length))"

    output *= "\n Mixed layer depth: $(summary(bgc.mixed_layer_depth))"

    output *= "\n Euphotic depth: $(summary(bgc.euphotic_depth))"

    output *= "\n Silicate climatology: $(summary(bgc.silicate_climatology))"

    output *= "\n Mixed layer mean diffusivity: $(summary(bgc.mean_mixed_layer_vertical_diffusivity))"

    output *= "\n Mixed layer mean light: $(summary(bgc.mean_mixed_layer_light))"

    output *= "\n Carbon chemistry: $(summary(bgc.carbon_chemistry))"

    output *= "\n Calcite saturation: $(summary(bgc.calcite_saturation))"

    output *= "\n Sinking velocities:"
    output *= "\n  Small particles: $(summary(bgc.sinking_velocities.POC))"
    output *= "\n  Large particles: $(summary(bgc.sinking_velocities.GOC))"

    print(io, output) 

    return nothing
end

summary(phyto::NanoAndDiatoms{<:Any, <:Any, FT}) where FT =
    string("NanoAndDiatoms{$(summary(phyto.nano)), $(summary(phyto.diatoms)), $FT} - $(required_biogeochemical_tracers(phyto))")

summary(phyto::NanoAndDiatoms{<:MixedMondo, <:MixedMondo, FT}) where FT =
    string("MixedMondo-NanoAndDiatoms{$FT} - $(required_biogeochemical_tracers(phyto))")

summary(::MixedMondo) = string("MixedMondo")

summary(zoo::MicroAndMeso{<:QualityDependantZooplankton, <:QualityDependantZooplankton, FT}) where FT =
    string("QualityDependantZooplankton-MicroAndMeso{$FT} - $(required_biogeochemical_tracers(zoo))")

summary(zoo::MicroAndMeso{<:Any, <:Any, FT}) where FT =
    string("MicroAndMeso {$(summary(zoo.micro)), $(summary(zoo.meso)), $FT} - $(required_biogeochemical_tracers(zoo))")

summary(::QualityDependantZooplankton) = string("QualityDependantZooplankton")

summary(::DissolvedOrganicCarbon{FT}) where FT = string("DissolvedOrganicCarbon{$FT} - DOC")

summary(pom::TwoCompartementCarbonIronParticles{FT}) where FT = 
    string("TwoCompartementCarbonIronParticles{$FT} - $(required_biogeochemical_tracers(pom))")

summary(::NitrateAmmonia{FT}) where FT = string("NitrateAmmonia{$FT} - (NO₃, NH₄)")
summary(::SimpleIron{FT}) where FT = string("SimpleIron{$FT} - Fe")
summary(::Oxygen{FT}) where FT = string("Oxygen{$FT} - O₂ (mmol O₂ / m³)")
summary(::Silicate) = string("Silicate - Si")
summary(::Phosphate) = string("Phosphate - PO₄")
summary(::InorganicCarbon) = string("InorganicCarbon - (DIC, Alk)")

summary(::ModelLatitude) = string("ModelLatitude")
summary(lat::PrescribedLatitude{FT}) where FT = string("PrescribedLatitude{FT} $(lat.latitude)°")

# TODO: add show methods
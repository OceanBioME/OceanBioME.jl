summary(::PISCES) = string("PISCES biogeochemical model") 

function show(io::IO, bgc::PISCES)

    FT = typeof(bgc.background_shear)   

    output = "Pelagic Interactions Scheme for Carbon and Ecosystem Studies (PISCES) model {$FT}"

    output *= "\n Nanophytoplankton: $(summary(bgc.nanophytoplankton))"

    output *= "\n Diatoms: $(summary(bgc.diatoms))"

    output *= "\n Microzooplankton: $(summary(bgc.microzooplankton))"

    output *= "\n Mesozooplankton: $(summary(bgc.mesozooplankton))"

    output *= "\n Microzooplankton: $(summary(bgc.microzooplankton))"

    output *= "\n Dissolved organic matter: $(summary(bgc.dissolved_organic_matter))"

    output *= "\n Particulate organic matter: $(summary(bgc.particulate_organic_matter))"

    output *= "\n Nitrogen: $(summary(bgc.nitrogen))"

    output *= "\n Iron: $(summary(bgc.iron))"

    output *= "\n Silicate: $(summary(bgc.silicate))"

    output *= "\n Oxygen: $(summary(bgc.oxygen))"

    output *= "\n Phosphate: $(summary(bgc.phosphate))"

    output *= "\n Calcite: $(summary(bgc.calcite))"

    output *= "\n Carbon system: $(summary(bgc.carbon_system))"

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

function summary(phyto::MixedMondoPhytoplankton{<:Any, <:Any, FT}) where FT
    growth_rate = summary(phyto.growth_rate)
    nutrient_limitation = summary(phyto.nutrient_limitation)

    names = "\n  I (mmol C / m³), IChl (mg Chl / m³), IFe (μmol Fe / m³)"

    if phyto.nutrient_limitation.silicate_limited
        names *= ", ISi (mmol Si / m³)"
    end
    
    return string("MixedMondoPhytoplankton{$(growth_rate), $(nutrient_limitation), $FT} - "*names)
end

summary(::NutrientLimitedProduction) = string("NutrientLimitedProduction")
summary(::GrowthRespirationLimitedProduction) = string("NutrientLimitedProduction")

summary(::NitrogenIronPhosphateSilicateLimitation) = string("NitrogenIronPhosphateSilicateLimitation")

summary(::Zooplankton{FT}) where FT = string("Zooplankton{$FT} - I (mmol C / m³)")

summary(::DissolvedOrganicMatter{FT}) where FT = string("DissolvedOrganicMatter{$FT} - DOC (mmol C / m³)")
summary(::TwoCompartementParticulateOrganicMatter{FT}) where FT = 
    string("TwoCompartementParticulateOrganicMatter{$FT} - 
  POC (mmol C / m³), GOC (mmol C / m³), SFe (μmol Fe / m³), BFe (μmol Fe / m³), PSi (mmol Si / m³)")

summary(::NitrateAmmonia{FT}) where FT = string("NitrateAmmonia{$FT} - NO₃ (mmol N / m³), NH₄(mmol N / m³)")
summary(::SimpleIron{FT}) where FT = string("SimpleIron{$FT} - Fe (μmol Fe / m³)")
summary(::Oxygen{FT}) where FT = string("Oxygen{$FT} - O₂ (mmol O₂ / m³)")
summary(::Silicate) = string("Silicate - Si (mmol Si / m³)")
summary(::Phosphate) = string("Phosphate - PO₄ (mmol P / m³)")
summary(::Calcite{FT}) where FT = string("Calcite{$FT} - CaCO₃ (mmol C / m³)")
summary(::CarbonateSystem) = string("CarbonateSystem - DIC (mmol C / m³), Alk (mequiv / m³)")

summary(::ModelLatitude) = string("ModelLatitude")
summary(lat::PrescribedLatitude{FT}) where FT = string("PrescribedLatitude $(lat.latitude)° {FT}")

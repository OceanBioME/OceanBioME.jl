function MixedMondoNanoAndDiatoms(FT = Float64; 
                                  nano = MixedMondo(FT; growth_rate = GrowthRespirationLimitedProduction{FT}(dark_tollerance = 3days),
                                                    nutrient_limitation = 
                                                        NitrogenIronPhosphateSilicateLimitation{FT}(minimum_ammonium_half_saturation = 0.013,
                                                                                                    minimum_nitrate_half_saturation = 0.13, 
                                                                                                    minimum_phosphate_half_saturation = 0.8,
                                                                                                    silicate_limited = false),
                                                    blue_light_absorption = convert(FT, 2.1), 
                                                    green_light_absorption = convert(FT, 0.42), 
                                                    red_light_absorption = convert(FT, 0.4),
                                                    maximum_quadratic_mortality = convert(FT, 0.0),
                                                    maximum_chlorophyll_ratio = convert(FT, 0.033),
                                                    half_saturation_for_iron_uptake = convert(FT, 1.0)),
                                  diatoms = MixedMondo(FT; growth_rate = GrowthRespirationLimitedProduction{FT}(dark_tollerance = 4days),
                                                       nutrient_limitation = 
                                                           NitrogenIronPhosphateSilicateLimitation{FT}(minimum_ammonium_half_saturation = 0.039,
                                                                                                       minimum_nitrate_half_saturation = 0.39, 
                                                                                                       minimum_phosphate_half_saturation = 2.4,
                                                                                                       silicate_limited = true),
                                                       blue_light_absorption = convert(FT, 1.6), 
                                                       green_light_absorption = convert(FT, 0.69), 
                                                       red_light_absorption = convert(FT, 0.7),
                                                       maximum_quadratic_mortality = convert(FT, 0.03/day),
                                                       maximum_chlorophyll_ratio = convert(FT, 0.05),
                                                       half_saturation_for_iron_uptake = convert(FT, 3.0)))
 
    return NanoAndDiatoms(; nano, diatoms)
end

const NANO_PHYTO = Union{Val{:P}, Val{:PChl}, Val{:PFe}}
const DIATOMS    = Union{Val{:D}, Val{:DChl}, Val{:DFe}, Val{:DSi}}

@inline phytoplankton_concentrations(::NANO_PHYTO, i, j, k, fields) = @inbounds fields.P[i, j, k], fields.PChl[i, j, k], fields.PFe[i, j, k]
@inline phytoplankton_concentrations(::DIATOMS,    i, j, k, fields) = @inbounds fields.D[i, j, k], fields.DChl[i, j, k], fields.DFe[i, j, k]

@inline carbon_name(::NANO_PHYTO) = Val(:P)
@inline carbon_name(::DIATOMS) = Val(:D)

@inline parameterisation(::NANO_PHYTO, phyto::NanoAndDiatoms) = phyto.nano
@inline parameterisation(::DIATOMS,    phyto::NanoAndDiatoms) = phyto.diatoms

# I think these could be abstracted more so that we have a few simple functions in nano_and_diatoms
# and most only exist in  `mixed_mondo.jl`
# also maybe should be dispatched on PISCES{NanoAndDiatoms{MixedMondo, MixedMondo}}
@inline function (bgc::PISCES{<:NanoAndDiatoms})(i, j, k, grid, val_name::Union{Val{:P}, Val{:D}}, clock, fields, auxiliary_fields)
    phyto = parameterisation(val_name, bgc.phytoplankton)
    
    growth = carbon_growth(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    linear_mortality, quadratic_mortality = mortality(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    death = (linear_mortality + quadratic_mortality)

    getting_grazed = grazing(bgc.zooplankton, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return growth - death - getting_grazed
end

@inline function (bgc::PISCES{<:NanoAndDiatoms})(i, j, k, grid, val_name::Union{Val{:PChl}, Val{:DChl}}, clock, fields, auxiliary_fields)
    phyto = parameterisation(val_name, bgc.phytoplankton)
    
    growth = chlorophyll_growth(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    I, IChl, IFe = phytoplankton_concentrations(val_name, i, j, k, fields)

    θChl = IChl / (12 * I + eps(0.0))

    linear_mortality, quadratic_mortality = mortality(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    death = (linear_mortality + quadratic_mortality)

    getting_grazed = grazing(bgc.zooplankton, carbon_name(val_name), i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    return growth - (death + getting_grazed) * θChl * 12
end

@inline function (bgc::PISCES{<:NanoAndDiatoms})(i, j, k, grid, val_name::Union{Val{:PFe}, Val{:DFe}}, clock, fields, auxiliary_fields)
    phyto = parameterisation(val_name, bgc.phytoplankton)
    
    growth = iron_growth(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    I, IChl, IFe = phytoplankton_concentrations(val_name, i, j, k, fields)

    θFe = IFe / (I + eps(0.0))

    linear_mortality, quadratic_mortality = mortality(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    death = (linear_mortality + quadratic_mortality)

    getting_grazed = grazing(bgc.zooplankton, carbon_name(val_name), i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    return growth - (death + getting_grazed) * θFe
end

@inline function (bgc::PISCES{<:NanoAndDiatoms})(i, j, k, grid, val_name::Val{:DSi}, clock, fields, auxiliary_fields)
    phyto = parameterisation(val_name, bgc.phytoplankton)
    
    growth = silicate_growth(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    D = @inbounds fields.D[i, j, k]
    DSi = @inbounds fields.DSi[i, j, k]

    θSi = DSi / (D + eps(0.0))

    linear_mortality, quadratic_mortality = mortality(phyto, val_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    death = (linear_mortality + quadratic_mortality)

    getting_grazed = grazing(bgc.zooplankton, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields) 

    return growth - (death + getting_grazed) * θSi
end
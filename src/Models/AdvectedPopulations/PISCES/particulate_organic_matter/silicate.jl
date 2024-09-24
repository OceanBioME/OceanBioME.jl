@inline function (bgc::TwoCompartementPOCPISCES)(i, j, k, grid, val_name::Val{:PSi}, clock, fields, auxiliary_fields)
    # this generalisation still assumes that we're only getting PSi from phytoplankton being grazed, will need changes if zooplankton get Si compartment

    phytoplankton_production = particulate_silicate_production(bgc.phytoplankton, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    dissolution = particulate_silicate_dissolution(bgc.particulate_organic_matter, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    return phytoplankton_production - dissolution
end

@inline function particulate_silicate_dissolution(poc::TwoCompartementCarbonIronParticles, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    PSi = @inbounds fields.PSi[i, j, k]
    Si  = @inbounds  fields.Si[i, j, k]

    T = @inbounds fields.T[i, j, k]
    
    λₗ = poc.fast_dissolution_rate_of_silicate
    λᵣ = poc.slow_dissolution_rate_of_silicate

    χ = particulate_silicate_liable_fraction(poc, i, j, k, grid, bgc, clock, fields, auxiliary_fields)

    λ₀ = χ * λₗ + (1 - χ) * λᵣ

    equilibrium_silicate = 10^(6.44 - 968 / (T + 273.15))
    silicate_saturation  = (equilibrium_silicate - Si) / equilibrium_silicate

    λ = λ₀ * (0.225 * (1 + T/15) * silicate_saturation + 0.775 * ((1 + T/400)^4 * silicate_saturation)^9)

    return λ * PSi # assuming the Diss_Si is typo in Aumont 2015, consistent with Aumont 2005
end

@inline function particulate_silicate_liable_fraction(poc, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    χ₀ = poc.base_liable_silicate_fraction
    λₗ = poc.fast_dissolution_rate_of_silicate
    λᵣ = poc.slow_dissolution_rate_of_silicate

    z = znode(i, j, k, grid, Center(), Center(), Center())

    zₘₓₗ = @inbounds auxiliary_fields.zₘₓₗ[i, j, k]
    zₑᵤ  = @inbounds  auxiliary_fields.zₑᵤ[i, j, k]

    # this isn't actually correct since wGOC isn't constant but nm
    wGOC = ℑzᵃᵃᶜ(i, j, k, grid, auxiliary_fields.wGOC) 

    zₘ = min(zₘₓₗ, zₑᵤ)

    return χ₀ * ifelse(z >= zₘ, 1, exp((λₗ - λᵣ) * (zₘ - z) / wGOC))
end

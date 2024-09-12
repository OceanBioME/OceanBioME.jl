@inline function (poc::TwoCompartementParticulateOrganicMatter)(bgc, ::Val{:PSi}, 
                                                                x, y, z, t,
                                                                P, D, Z, M, 
                                                                PChl, DChl, PFe, DFe, DSi,
                                                                DOC, POC, GOC, 
                                                                SFe, BFe, PSi, 
                                                                NO₃, NH₄, PO₄, Fe, Si, 
                                                                CaCO₃, DIC, Alk, 
                                                                O₂, T, 
                                                                zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)

    # diatom grazing
    _, _, microzooplankton_grazing = specific_grazing(bgc.microzooplankton, P, D, Z, POC)
    _, _, mesozooplankton_grazing  = specific_grazing(bgc.mesozooplankton , P, D, Z, POC)

    diatom_grazing = (microzooplankton_grazing * Z + mesozooplankton_grazing * M) * DSi / D

    # diatom mortality
    diatom_linear_mortality, diatom_quadratic_mortality = mortality(bgc.diatoms, bgc, z, D, zₘₓₗ)

    diatom_mortality = (diatom_linear_mortality + diatom_quadratic_mortality) * DSi / D

    # dissolution
    dissolution = particulate_silicate_dissolution(poc, bgc, x, y, z, PSi, Si, T, zₘₓₗ, zₑᵤ)

    return (diatom_grazing + diatom_mortality - dissolution)
end

@inline function particulate_silicate_dissolution(poc, bgc, x, y, z, PSi, Si, T, zₘₓₗ, zₑᵤ)
    λₗ = poc.fast_dissolution_rate_of_silicate
    λᵣ = poc.slow_dissolution_rate_of_silicate

    χ = particulate_silicate_liable_fraction(poc, bgc, x, y, z, zₘₓₗ, zₑᵤ)

    λ₀ = χ * λₗ + (1 - χ) * λᵣ

    equilibrium_silicate = 10^(6.44 - 968 / (T + 273.15))
    silicate_saturation  = (equilibrium_silicate - Si) / Si

    λ = λ₀ * (0.225 * (1 + T/15) * silicate_saturation + 0.775 * ((1 + T/400)^4 * silicate_saturation)^9)

    return λ * PSi # assuming the Diss_Si is typo in Aumont 2015, consistent with Aumont 2005
end

@inline function particulate_silicate_liable_fraction(poc, bgc, x, y, z, zₘₓₗ, zₑᵤ)
    χ₀ = poc.base_liable_silicate_fraction
    λₗ = poc.fast_dissolution_rate_of_silicate
    λᵣ = poc.slow_dissolution_rate_of_silicate
    grid = bgc.sinking_velocities.grid

    w = particle_sinking_speed(x, y, z, grid,  bgc.sinking_velocities.GOC.w)

    zₘ = min(zₘₓₗ, zₑᵤ)

    return χ₀ * ifelse(z >= zₘ, 1, exp((λₗ - λᵣ) * (zₘ - z) / w))
end

@kwdef struct Zooplankton{FT}
    temperature_sensetivity :: FT = 1.079
    maximum_grazing_rate :: FT

    preference_for_nanophytoplankton :: FT
    preference_for_diatoms :: FT
    preference_for_particulates :: FT
    preference_for_zooplankton :: FT

    food_threshold_concentration :: FT = 0.3
    specific_food_thresehold_concentration :: FT = 0.001

    grazing_half_saturation :: FT = 20.0

    maximum_flux_feeding_rate :: FT = 2.0e-3 # is this correct, paper says 2e3

    iron_ratio :: FT = 10^-3 # units?

    maximum_growth_efficiency :: FT
    non_assililated_fraction :: FT = 0.3

    mortality_half_saturation :: FT = 0.2
    quadration_mortality :: FT
    linear_mortality :: FT

    dissolved_excretion_fraction :: FT = 0.6
end

@inline zooplankton_concentration(::Val{Z}, Z, M) = Z
@inline zooplankton_concentration(::Val{M}, Z, M) = M

@inline function specific_grazing(zoo::Zooplankton, P, D, Z, POC)
    g₀   = zoo.maximum_grazing_rate
    b    = zoo.temperature_sensetivity
    pP   = zoo.preference_for_nanophytoplankton
    pD   = zoo.preference_for_diatoms
    pPOC = zoo.preference_for_particulates
    pZ   = zoo.preference_for_zooplankton
    J    = zoo.specific_food_thresehold_concentration
    K    = zoo.grazing_half_saturation

    food_threshold_concentration = zoo.food_threshold_concentration

    base_grazing_rate = g₀ * b ^ T

    food_availability = pP * P + pD * D + pPOC * POC + pZ * Z

    weighted_food_availability = pP * max(0, P - J) + pD * max(0, D - J) + pPOC * max(0, POC - J) + pZ * max(0, Z - J)

    concentration_limited_grazing = max(0, weighted_food_availability - min(weighted_food_availability / 2, food_threshold_concentration))

    total_specific_grazing = base_grazing_rate * concentration_limited_grazing / (K + food_availability) 

    phytoplankton_grazing = pP * max(0, P - J)     * total_specific_grazing / weighted_food_availability
    diatom_grazing        = pD * max(0, D - J)     * total_specific_grazing / weighted_food_availability
    particulate_grazing   = pPOC * max(0, POC - J) * total_specific_grazing / weighted_food_availability
    zooplankton_grazing   = pZ * max(0, Z - J)     * total_specific_grazing / weighted_food_availability

    return total_specific_grazing, phytoplankton_grazing, diatom_grazing, particulate_grazing, zooplankton_grazing
end

@inline function specific_flux_feeding(zoo::Zooplankton, POC, w_field, grid)
    g₀ = zoo.maximum_flux_feeding_rate
    b  = zoo.temperature_sensetivity

    base_flux_feeding_rate = g₀ * b ^ T

    # hopeflly this works on GPU
    w = particle_sinking_speed(x, y, z, grid, w_field)

    return base_flux_feeding_rate * w * POC
end

@inline function (zoo::Zooplankton)(bgc, val_name::Union{Val{:Z}, Val{:M}}, 
                                    x, y, z, t,
                                    P, D, Z, M, 
                                    PChl, DChl, PFe, DFe, DSi, 
                                    DOC, POC, GOC, 
                                    SFe, BFe, PSi, 
                                    NO₃, NH₄, PO₄, Fe, Si, 
                                    CaCO₃, DIC, Alk, 
                                    O₂, T, 
                                    zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, PAR, PAR₁, PAR₂, PAR₃)

    I = zooplankton_concentration(val_name, Z, M)

    # grazing
    total_specific_grazing, gP, gD, gPOC, gZ = specific_grazing(zoo, P, D, Z, POC)

    grazing = total_specific_grazing * I

    # growth efficiency - don't need the nitrogen term since N/C is the same for everything
    θFe = zoo.iron_ratio
    e₀  = zoo.maximum_growth_efficiency
    σ   = zoo.non_assililated_fraction

    iron_grazing = PFe / P * gP + DFe / D * gD + SFe / POC * gPOC + θFe * gZ

    iron_grazing_ratio = iron_grazing / (θFe * total_specific_grazing)

    food_quality = min(1, iron_grazing_ratio)

    growth_efficiency = food_quality * min(e₀, (1 - σ) * iron_grazing_ratio)

    # flux feeding
    grid = bgc.sinking_velocities.grid
    small_flux_feeding = specific_flux_feeding(zoo, POC, bgc.sinking_velocities.POC.w, grid)
    large_flux_feeding = specific_flux_feeding(zoo, GOC, bgc.sinking_velocities.GOC.w, grid)

    flux_feeding = (small_flux_feeding + large_flux_feeding) * I

    # grazing mortality
    specific_grazing_mortality = grazing_mortality(val_name, bgc.mesozooplankton, P, D, Z, POC)

    grazing_mortality = specific_grazing_mortality * M

    # mortality
    b  = zoo.temperature_sensetivity
    m₀ = zoo.quadratic_mortality
    Kₘ = zoo.mortality_half_saturation
    r  = zoo.linear_mortality

    temperature_factor = b^T

    concentration_factor = Z / (Z + Kₘ)

    mortality = temperature_factor * Z * (m₀ * Z + r * (concentration_factor + 3 * anoxia_factor(bgc, O₂)))



    return growth_efficiency * (grazing + flux_feeding) - grazing_mortality - mortality
end

@inline function nanophytoplankton_grazing(zoo::Zooplankton, P, D, Z, POC) 
    _, g = specific_grazing(zoo, P, D, Z, POC)

    return g
end

@inline function diatom_grazing(zoo::Zooplankton, P, D, Z, POC) 
    _, _, g = specific_grazing(zoo, P, D, Z, POC)

    return g
end

@inline function particulate_grazing(zoo::Zooplankton, P, D, Z, POC) 
    _, _, _, g = specific_grazing(zoo, P, D, Z, POC)

    return g
end

@inline function zooplankton_grazing(zoo::Zooplankton, P, D, Z, POC) 
    _, _, _, _, g = specific_grazing(zoo, P, D, Z, POC)

    return g
end

@inline specific_small_flux_feeding(zoo::Zooplankton, bgc, POC, GOC) =
    specific_flux_feeding(zoo, POC, bgc.sinking_velocities.POC.w, grid)

@inline specific_large_flux_feeding(zoo::Zooplankton, bgc, POC, GOC) =
    specific_flux_feeding(zoo, POC, bgc.sinking_velocities.POC.w, grid)

@inline grazing_mortality(val_name, zoo, P, D, Z, POC) = 0
@inline grazing_mortality(::Val{:Z}, zoo, P, D, Z, POC) = zooplankton_grazing(zoo, P, D, Z, POC)

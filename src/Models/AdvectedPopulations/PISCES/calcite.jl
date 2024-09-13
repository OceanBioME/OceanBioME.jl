@kwdef struct Calcite{FT}
          base_rain_ratio :: FT = 0.3
    base_dissolution_rate :: FT = 0.197/day
     dissolution_exponent :: FT = 1.0
end

@inline function (calcite::Calcite)(bgc, ::Val{:CaCO₃}, 
                                    x, y, z, t,
                                    P, D, Z, M, 
                                    PChl, DChl, PFe, DFe, DSi,
                                    DOC, POC, GOC, 
                                    SFe, BFe, PSi, 
                                    NO₃, NH₄, PO₄, Fe, Si, 
                                    CaCO₃, DIC, Alk, 
                                    O₂, T, 
                                    zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    production = calcite_production(calacite, bgc, P, D, PChl, PFe, Z, M, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    dissolution = calcite_dissolution(calcite, CaCO₃, Ω)
    
    return production - dissolution
end

@inline function calcite_production(calcite, bgc, P, D, PChl, PFe, Z, M, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)
    R = rain_ratio(calcite, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)

    microzooplankton = specific_calcite_grazing_loss(bgc.microzooplankton, P, D, Z, POC) * Z
    mesozooplankton  = specific_calcite_grazing_loss(bgc.mesozooplankton, P, D, Z, POC) * M

    linear_mortality, quadratic_mortality = mortality(bgc.nanophytoplankton, bgc, z, I, zₘₓₗ)

    mortality = 0.5 * (linear_mortality + quadratic_mortality)

    return R * (microzooplankton + mesozooplankton + mortality)
end

# should this be in the particles thing?
@inline function rain_ratio(calcite::Calcite, bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, PAR)
    r = calcite.base_rain_ratio

    # assuming this is a type in Aumont 2015 based on Aumont 2005
    L, = bgc.nanophytoplankton.nutrient_limitation(bgc, P, PChl, PFe, NO₃, NH₄, PO₄, Fe, Si, Si′)

    L_CaCO₃ = L # maybe this is wrong, can't find the reference, others had it as min(Lₙᴾ, concentration_limitation(Fe, 6e-11), concentration_limitation(PO₄, Kₙₕ₄ᴾ))

    phytoplankton_concentration_factor = max(1, P / 2)

    low_light_factor = max(0, PAR - 1) / (4 + PAR)
    high_light_factor = 30 / (30 + PAR)

    low_temperature_factor = max(0, T / (T + 0.1)) # modified from origional as it goes negative and does not achieve goal otherwise
    high_temperature_factor = 1 + exp(-(T - 10)^2 / 25)

    depth_factor = min(1, -50/zₘₓₗ)

    return r * L_CaCO₃ * phytoplankton_concentration_factor * low_light_factor * high_light_factor * low_temperature_factor * high_temperature_factor * depth_factor
end

@inline function calcite_dissolution(calcite, CaCO₃, Ω)
    λ   = calcite.base_dissolution_rate
    nca = calcite.dissolution_exponent

    ΔCaCO₃ = max(0, 1 - Ω)

    return λ * ΔCaCO₃ ^ nca * CaCO₃
end

# there is a typo in https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2010JC006446 so I am not sure the first term is correct, but this makes sense
@inline function (bgc::LOBSTER)(::Val{:OXY}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR)
    γ = bgc.phytoplankton_exudation_fraction
    μₚ = bgc.maximum_phytoplankton_growthrate
    kₚₐᵣ = bgc.light_half_saturation
    ψ = bgc.nitrate_ammonia_inhibition
    kₙₒ₃ = bgc.nitrate_half_saturation
    ROᵖ = bgc.respiraiton_oxygen_nitrogen_ratio
    ROⁿ = bgc.nitrifcation_oxygen_nitrogen_ratio
    μₙ = bgc.nitrifcaiton_rate

    return (μₚ * Lₚₐᵣ(PAR, kₚₐᵣ) * Lₙₒ₃(NO₃, NH₄, ψ, kₙₒ₃) * ROᵖ * P
            - (ROᴾ - ROⁿ)*bgc(Val(:NH₄), x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
            - ROᴾ * μₙ * NH₄)
end
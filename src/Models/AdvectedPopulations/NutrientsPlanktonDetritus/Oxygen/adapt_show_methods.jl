function Adapt.adapt_structure(to, ox::RedoxOxygen{FT}) where FT
    oxygen_consumption_scale_factor = adapt(to, ox.oxygen_consumption_scale_factor)

    return RedoxOxygen{FT, typeof(oxygen_consumption_scale_factor)}(
        ox.nitrate_carbon_to_oxygen,
        ox.remineralisation_carbon_to_oxygen,
        ox.diazotroph_carbon_to_oxygen,
        ox.denitrification_carbon_to_nitrogen,
        ox.nitrification_rate,
        ox.nitrification_par_limit,
        ox.oxygen_minimum,
        ox.oxygen_minimum_range,
        ox.fastest_depletion_timescale,
        oxygen_consumption_scale_factor
    )
end

Base.summary(::RedoxOxygen) =
    "RedoxOxygen (:O₂; nitrification NH₄→NO₃, denitrification NO₃→N₂, N-source-dependent O₂ production)"

function Base.show(io::IO, ox::RedoxOxygen)
    msg  = summary(ox) * "\n"
    msg *= "├── C:O₂ (nitrate/remineralisation/diazotroph): "
    msg *= "$(round(ox.nitrate_carbon_to_oxygen, sigdigits = 4))/"
    msg *= "$(round(ox.remineralisation_carbon_to_oxygen, sigdigits = 4))/"
    msg *= "$(round(ox.diazotroph_carbon_to_oxygen, sigdigits = 4))\n"
    msg *= "├── nitrification: $(ox.nitrification_rate)/s below $(ox.nitrification_par_limit) W/m²\n"
    msg *= "└── oxygen minimum: $(ox.oxygen_minimum) ± $(ox.oxygen_minimum_range) mmol/m³"

    print(io, msg)

    return nothing
end

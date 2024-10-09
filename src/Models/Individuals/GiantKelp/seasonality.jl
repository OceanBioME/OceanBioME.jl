using Oceananigans.Units

using OceanBioME: CBMDayLength # from utils

struct NotSeasonal end

(::NotSeasonal)(args...) = 1

@kwdef struct RateOfDayLengthChange{FT, DL, PC}
           photoperiod_1 :: FT = 0.85 # stupid names, no idea what they actually mean
           photoperiod_2 :: FT = 0.3
                latitude :: FT # TODO: make this a function, also somehow get the latitude into the particles
              day_length :: DL = CBMDayLength()
  peak_day_length_change :: PC = find_maximum_day_length_change(day_length, latitude) # don't think theres an analytical solution for this
end

function find_maximum_day_length_change(day_length, φ)
    day_lengths = map(t->day_length(t, φ), [0:1day:365day;])

    ΔL = diff(day_lengths)

    return findmax(abs, ΔL)[1]
end

# the kelp magically know what day of the year it is!!!
@inline function (seasonality::RateOfDayLengthChange)(t)
    φ = seasonality.latitude
    a₁ = seasonality.photoperiod_1
    a₂ = seasonality.photoperiod_2
    ΔLₚ = seasonality.peak_day_length_change # if they're going to move this will have to vary with latitude

    Lₙ = seasonality.day_length(t, φ)
    Lₙ₋₁ = seasonality.day_length((floor(Int, mod(t, 365days)/day) - 1)*day, φ)

    λ = (Lₙ - Lₙ₋₁) / ΔLₚ

    return a₁ * (1 + sign(λ) * √abs(λ)) + a₂
end

using Oceananigans.Units

include("callablevalue.jl")
include("timestep.jl")
include("negative_tracers.jl")
include("sinking_velocity_fields.jl")
include("solvers.jl")
include("unwrapvaluefields.jl")

"""
    (day_length::CBMDayLength)(t, φ)

Returns the length of day in seconds at the latitude `φ`, `t` seconds after the start of the year.
"""
@kwdef struct CBMDayLength{FT}
    day_length_coefficient :: FT = 0.833
end

# TODO: add methods for DateTime times etc
@inline function (day_length::CBMDayLength{FT})(t, φ) where FT
    # as per Forsythe et al., 1995 (https://doi.org/10.1016/0304-3800(94)00034-F)
    p = day_length.day_length_coefficient

    # day of year
    J = floor(Int, mod(t, 365days)/day)

    # revolution angle
    θ = 0.216310 + 2 * atan(0.9671396 * tan(0.00860 * (J - 186)))

    # solar declination
    ϕ = asind(0.39795 * cos(θ))

    L = max(-one(t), min(one(t), (sind(p) + sind(φ)*sind(ϕ)) / (cosd(φ) * cosd(ϕ))))

    return convert(FT, (24 - 24 / 180 * acosd(L)) * hour)
end

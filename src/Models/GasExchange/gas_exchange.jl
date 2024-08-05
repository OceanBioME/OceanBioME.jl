"""
    GasExchange

`GasExchange` returns the air-sea flux of a gas betwen `water_concentration` and
`air_concentration` with a `transfer_velocity` computed from the temperature 
(provided later), and the `wind_speed`.

### TODO: UPDATE
`transfer_velocity` should behave as a function of wind speed and temperature (i.e.
`k(u, T)`), `water_concentration` a function of `c(x, y, t, T, field_dependencies...)`,
and `air_concentration` of `ac(x, y, t)`, where `field_dependencies` are specified.
"""
struct GasExchange{WS, TV, WC, AC} <: Function
             wind_speed :: WS
      transfer_velocity :: TV
    water_concentration :: WC
      air_concentration :: AC
end

@inline function (g::GasExchange)(i, j, grid, clock, model_fields)
    T = @inbounds model_fields.T[i, j, grid.Nx]

    u₁₀ = surface_value(g.wind_speed, i, j, grid, clock)

    k = g.transfer_velocity(u₁₀, T)

    air_concentration = surface_value(g.air_concentration, i, j, grid, clock)

    water_concentration = surface_value(g.water_concentration, i, j, grid, clock, model_fields)

    return k * (water_concentration - air_concentration)
end

Adapt.adapt_structure(to, g::GasExchange) = GasExchange(adapt(to, g.wind_speed),
                                                        adapt(to, g.transfer_velocity),
                                                        adapt(to, g.water_concentration),
                                                        adapt(to, g.air_concentration),
                                                        nothing)

summary(::GasExchange{WS, TV, WC}) where {WS, TV, WC} = "Air-sea `GasExchange` model for $(nameof(WC))"
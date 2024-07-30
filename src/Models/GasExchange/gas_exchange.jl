"""
    GasExchange

`GasExchange` returns the air-sea flux of a gas betwen `water_concentration` and
`air_concentration` with a `transfer_velocity` computed from the temperature 
(provided later), and the `wind_speed`.

`transfer_velocity` should behave as a function of wind speed and temperature (i.e.
`k(u, T)`), `water_concentration` a function of `c(x, y, t, T, field_dependencies...)`,
and `air_concentration` of `ac(x, y, t)`, where `field_dependencies` are specified.
"""
struct GasExchange{WS, TV, WC, AC, FD} <: Function
             wind_speed :: WS
      transfer_velocity :: TV
    water_concentration :: WC
      air_concentration :: AC
     field_dependencies :: FD
end

@inline function (g::GasExchange)(x, y, t, T, args...)
    k = g.transfer_velocity(surface_value(g.wind_speed, x, y, t), T)

    conc = surface_value(g.water_concentration, x, y, t, T, args...)

    air_conc = surface_value(g.air_concentration, x, y, t)

    return k * (conc - air_conc)
end

Adapt.adapt_structure(to, g::GasExchange) = GasExchange(adapt(to, g.wind_speed),
                                                        adapt(to, g.transfer_velocity),
                                                        adapt(to, g.water_concentration),
                                                        adapt(to, g.air_concentration),
                                                        nothing)

summary(::GasExchange{WS, TV, WC}) where {WS, TV, WC} = "Air-sea gas exchange model for $(nameof(WC))"
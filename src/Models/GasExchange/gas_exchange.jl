struct GasExchange{WS, TV, WC, AC, FD, FM} <: Function
             wind_speed :: WS
      transfer_velocity :: TV
    water_concentration :: WC
      air_concentration :: AC
     field_dependencies :: FD
 field_dependencies_map :: FM
end

@inline function (g::GasExchange)(x, y, t, args...)
    k = g.transfer_velocity(surface_value(g.wind_speed, x, y, t), args[g.field_dependencies_map.T])

    conc = surface_value(g.water_concentration, g.field_dependencies_map, x, y, t, args...)

    air_conc = surface_value(g.air_concentration, x, y, t)

    return k * (conc - air_conc)
end

Adapt.adapt_structure(to, g::GasExchange) = GasExchange(adapt(to, g.wind_speed),
                                                        adapt(to, g.transfer_velocity),
                                                        adapt(to, g.water_concentration),
                                                        adapt(to, g.air_concentration),
                                                        adapt(to, g.field_dependencies),
                                                        adpat(to, g.field_dependencies_map))
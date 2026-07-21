import Oceananigans.BoundaryConditions: regularize_boundary_condition, DiscreteBoundaryFunction

carbon_dioxide_exchange_bc = DiscreteBoundaryFunction{<:Any, <:GasExchange{<:Any, <:Any, <:CarbonDioxideConcentration{}}}

function regularize_boundary_condition(wrapped_condition::carbon_dioxide_exchange_bc, 
                                       grid,
                                       loc, 
                                       dim,
                                       side,
                                       dic_name)

    condition = wrapped_condition.func

    @info dic_name
   
    if !is_materialised(condition.water_concentration)
        identifier = @show last(string(dic_name))

        if isnumeric(identifier)
            DIC = Symbol(:DIC, identifier)
            Alk = Symbol(:Alk, identifier)
        else
            DIC = :DIC
            Alk = :Alk
        end

        wc = condition.water_concentration

        water_concentration = 
            CarbonDioxideConcentration{DIC, Alk}(wc.carbon_chemistry,
                                                 wc.first_virial_coefficient,
                                                 wc.cross_virial_coefficient,
                                                 wc.air_pressure,
                                                 wc.silicate_and_phosphate_names)

        return GasExchange(condition.wind_speed,
                           condition.transfer_velocity,
                           water_concentration,
                           condition.air_concentration)
    else
        return condition
    end
end

is_materialised(::CarbonDioxideConcentration) = true
is_materialised(::CarbonDioxideConcentration{<:Nothing}) = false

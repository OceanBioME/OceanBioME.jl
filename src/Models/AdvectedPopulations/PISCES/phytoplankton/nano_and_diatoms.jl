@kwdef struct NanoAndDiatoms{N, D, FT}
               nano :: N
            diatoms :: D
    base_rain_ratio :: FT = 0.3
end

required_biogeochemical_tracers(phyto::NanoAndDiatoms) = (required_biogeochemical_tracers(phyto.nano, :P)...,
                                                          required_biogeochemical_tracers(phyto.diatoms, :D)...)

@inline dissolved_exudate(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (dissolved_exudate(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + dissolved_exudate(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline uptake(phyto::NanoAndDiatoms, val_uptake_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (uptake(phyto.nano, Val(:P), val_uptake_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + uptake(phyto.diatoms, Val(:D), val_uptake_name, i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline function nitrogen_availability_limitation(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields)
    _, _, _, LN = phyto.nano.nutrient_limitation(Val(:P), i, j, k, grid, bgc, phyto.nano, clock, fields, auxiliary_fields)

    return LN
end

@inline base_production_rate(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    @inbounds base_production_rate(phyto.nano.growth_rate, fields.T[i, j, k])

@inline silicate_uptake(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (silicate_uptake(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + silicate_uptake(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields))

@inline total_production(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields, auxiliary_fields) =
    (total_production(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields, auxiliary_fields)
     + total_production(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields, auxiliary_fields))
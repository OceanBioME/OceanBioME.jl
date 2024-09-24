struct NanoAndDiatoms{N, D}
       nano :: N
    diatoms :: D
end

required_biogeochemical_tracers(phyto::NanoAndDiatoms) = (required_biogeochemical_tracers(phyto.nano, :P)...,
                                                          required_biogeochemical_tracers(phyto.diatoms, :D)...)

@inline dissolved_exudate(phyto::NanoAndDiatoms, bgc, i, j, k, grid, bgc, clock, fields) =
    (dissolved_exudate(phyto.nano, Val(:P), bgc, i, j, k, grid, bgc, clock, fields)
     + dissolved_exudate(phyto.diatoms, Val(:D), bgc, i, j, k, grid, bgc, clock, fields))

@inline uptake(phyto::NanoAndDiatoms, val_uptake_name, i, j, k, grid, bgc, clock, fields) =
    (uptake(phyto.nano, Val(:P), val_uptake_name, bgc, i, j, k, grid, bgc, clock, fields)
     + uptake(phyto.diatoms, Val(:D), val_uptake_name, bgc, i, j, k, grid, bgc, clock, fields))

@inline function nitrogen_availability_limitation(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields)
    _, _, _, LN = phyto.nano.nutrient_limitation(Val(:P), phyto.nano, i, j, k, grid, bgc, clock, fields)

    return LN
end

@inline base_production_rate(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields) =
    @inbounds base_production_rate(phyto.nano.growth_rate, fields.T[i, j, k])

@inline silicate_uptake(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields) =
    (silicate_uptake(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields)
     + silicate_uptake(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields))

@inline total_production(phyto::NanoAndDiatoms, i, j, k, grid, bgc, clock, fields) =
    (total_production(phyto.nano, Val(:P), i, j, k, grid, bgc, clock, fields)
     + total_production(phyto.diatoms, Val(:D), i, j, k, grid, bgc, clock, fields))
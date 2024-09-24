@inline function dissolved_exudate(phyto::MixedMondo, val_name, bgc, i, j, k, grid, bgc, clock, fields)
    δ  = phyto.exudated_fracton

    μI = total_production(phyto, val_name, bgc, i, j, k, grid, bgc, clock, fields)

    return δ * μI
end

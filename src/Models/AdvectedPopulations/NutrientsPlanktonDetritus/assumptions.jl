# the default assumptions

function dissolved_waste end
function solid_waste end
function calcite_dissolution end
function inorganic_waste end
function nutrient_uptake end

@inline nitrogen_ratio(i, j, k, grid, plankton, bgc, fields) = nitrogen_ratio(plankton, bgc)
@inline nitrogen_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = one(FT)

@inline carbon_ratio(i, j, k, grid, plankton, bgc, fields) = carbon_ratio(plankton, bgc)
@inline carbon_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = convert(FT, 106/16)

@inline phosphate_ratio(i, j, k, grid, plankton, bgc, fields) = phosphate_ratio(plankton, bgc)
@inline phosphate_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = convert(FT, 1/16)

@inline iron_ratio(i, j, k, grid, plankton, bgc, fields) = iron_ratio(plankton, bgc)
@inline iron_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = convert(FT, 0.0032/16) 

@inline silicon_ratio(i, j, k, grid, plankton, bgc, fields) = silicon_ratio(plankton, bgc)
@inline silicon_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = zero(FT)

@inline calcite_rain_ratio(i, j, k, grid, plankton, bgc, fields) = calcite_rain_ratio(plankton, bgc)
@inline calcite_rain_ratio(plankton, ::NutrientsPlanktonDetritus{FT}) where FT = zero(FT)

for (element, symbol) in pairs((nitrogen = :N, phosphate = :PO₄, iron = :Fe, silicon = :Si, carbon = :C))
    inorganic_waste_name  = Symbol(:inorganic_, element, :_waste)
    solid_waste_name      = Symbol(:solid_,     element, :_waste)
    dissolved_waste_name  = Symbol(:dissolved_, element, :_waste)
    ratio_name            = Symbol(element, :_ratio)
    @eval begin
        @inline $inorganic_waste_name(i, j, k, grid, plankton_or_detritus, bgc, fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            inorganic_waste(i, j, k, grid, plankton_or_detritus, bgc, fields, auxiliary_fields)

        @inline $solid_waste_name(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            solid_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

        @inline $dissolved_waste_name(i, j, k, grid, plankton, bgc, fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            dissolved_waste(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)

        @inline nutrient_uptake(i, j, k, grid, ::Val{$(QuoteNode(symbol))}, plankton, bgc, fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            nutrient_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
    end
end

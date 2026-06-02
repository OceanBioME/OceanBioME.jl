# the default assumptions

@inline nitrogen_ratio(i, j, k, grid, plankton, ::NutrientsPlanktonDetritus{FT}, fields) where FT = 
    one(FT)

@inline carbon_ratio(i, j, k, grid, plankton, ::NutrientsPlanktonDetritus{FT}, fields) where FT = 
    convert(FT, 106/16)

@inline phosphate_ratio(i, j, k, grid, plankton, ::NutrientsPlanktonDetritus{FT}, fields) where FT = 
    convert(FT, 1/16)

@inline iron_ratio(i, j, k, grid, plankton, ::NutrientsPlanktonDetritus{FT}, fields) where FT = 
    convert(FT, 0.0032/16) 

@inline silicon_ratio(i, j, k, grid, plankton, ::NutrientsPlanktonDetritus{FT}, fields) where FT = 
    zero(FT)

@inline calcite_rain_ratio(i, j, k, grid, plankton, ::NutrientsPlanktonDetritus{FT}, fields) where FT = 
    zero(FT)

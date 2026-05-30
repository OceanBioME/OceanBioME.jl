# the default assumptions

@inline nitrogen_ratio(plankton, ::NutrientsPlanktonDetritus{FT}, i, j, k, fields) where FT = 
    one(FT)

@inline carbon_ratio(plankton, ::NutrientsPlanktonDetritus{FT}, i, j, k, fields) where FT = 
    convert(FT, 106/16)

@inline phosphate_ratio(plankton, ::NutrientsPlanktonDetritus{FT}, i, j, k, fields) where FT = 
    convert(FT, 1/16)

@inline iron_ratio(plankton, ::NutrientsPlanktonDetritus{FT}, i, j, k, fields) where FT = 
    convert(FT, 0.0032/16) 

@inline silicon_ratio(plankton, ::NutrientsPlanktonDetritus{FT}, i, j, k, fields) where FT = 
    zero(FT)

@inline calcite_rain_ratio(plankton, ::NutrientsPlanktonDetritus{FT}, i, j, k, fields) where FT = 
    zero(FT)

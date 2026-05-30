# defaults
# defines solid_X_waste to be X_ratio * inorganic_waste and dissolved...
function solid_waste end
function dissolved_waste end

for (element, symbol) in pairs((nitrogen = :N, phosphate = :PO₄, iron = :Fe, silicon = :Si))
    solid_waste_name = Symbol(:solid_, element, :_waste)
    dissolved_waste_name = Symbol(:dissolved_, element, :_waste)
    ratio_name = Symbol(element, :_ratio)
    @eval begin
        @inline $solid_waste_name(plankton,
                                       bgc, i, j, k,
                                       fields, auxiliary_fields) = 
            $ratio_name(bgc.plankton, bgc, i, j, k, fields) *
            solid_waste(plankton, 
                        bgc, i, j, k,
                        fields, auxiliary_fields)

        @inline $dissolved_waste_name(plankton,
                                       bgc, i, j, k,
                                       fields, auxiliary_fields) = 
            $ratio_name(bgc.plankton, bgc, i, j, k, fields) *
            dissolved_waste(plankton, 
                            bgc, i, j, k,
                            fields, auxiliary_fields)
    end
end

@inline grazing(plankton, ::NPD{FT}, i, j, k, val_name, fields, auxiliary_fields) where FT = zero(FT)
@inline calcite_precipitation(::NPD{FT}, i, j, k, fields, auxiliary_fields) where FT = zero(FT)
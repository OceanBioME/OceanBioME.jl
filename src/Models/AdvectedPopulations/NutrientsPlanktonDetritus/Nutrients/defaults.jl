# defaults
# defines inorganic_X_waste to be X_ratio * inorganic_waste
# and nutrient_uptake for X element to be X_ratio * nutrient_uptake 

function inorganic_waste end

for (element, symbol) in pairs((nitrogen = :N, phosphate = :PO₄, iron = :Fe, silicon = :Si))
    inorganic_waste_name = Symbol(:inorganic_, element, :_waste)
    ratio_name = Symbol(element, :_ratio)
    
    @eval begin
        @inline $inorganic_waste_name(i, j, k, grid, 
                                      plankton_or_detritus, bgc,
                                      fields, auxiliary_fields) = 
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            inorganic_waste(i, j, k, grid, 
                            plankton_or_detritus, bgc, 
                            fields, auxiliary_fields)

        @inline nutrient_uptake(i, j, k, grid, 
                                plankton, bgc, 
                                ::Val{$(QuoteNode(symbol))},
                                fields, auxiliary_fields) =
            $ratio_name(i, j, k, grid, bgc.plankton, bgc, fields) *
            nutrient_uptake(i, j, k, grid, plankton, bgc, fields, auxiliary_fields)
    end
end

import OceanBioME: conserved_tracers

# possible conservations are N/P/Fe/Si/C/O
# neglect Si for now
# when a nutrient class is nothing it is implicitly conserved so we will not include that in the conservations
@inline function conserved_tracers(bgc::NutrientsPlanktonDetritus)
    conserved_elements = available_nutrients(bgc.nutrients)

    if !isnothing(bgc.inorganic_carbon)
        conserved_elements = (conserved_elements..., :carbon)
    end

    if !isnothing(bgc.oxygen)
        conserved_elements = (conserved_elements..., :oxygen)
    end

    conserved_tracers = NamedTuple()

    for element in conserved_elements
        element_tracers = group_element_tracers(bgc.nutrients, bgc, Val(element))
        element_tracers = merge(element_tracers, group_element_tracers(bgc.plankton, bgc, Val(element)))
        element_tracers = merge(element_tracers, group_element_tracers(bgc.detritus, bgc, Val(element)))
        element_tracers = merge(element_tracers, group_element_tracers(bgc.inorganic_carbon, bgc, Val(element)))
        element_tracers = merge(element_tracers, group_element_tracers(bgc.oxygen, bgc, Val(element)))

        conserved_tracers = merge(conserved_tracers, NamedTuple{(element, )}((element_tracers, ))) 
    end

    return conserved_tracers
end

# tracer conservations

group_element_tracers(Nutrients::Nutrients, bgc, val_element) = NamedTuple()
group_element_tracers(nutrients::Nutrients, bgc, val_element::Val{:nitrogen}) = 
    group_element_tracers(nutrients.nitrogen, bgc, val_element)
group_element_tracers(nutrients::Nutrients, bgc, val_element::Val{:phosphate}) = 
    group_element_tracers(nutrients.phosphate, bgc, val_element)
group_element_tracers(nutrients::Nutrients, bgc, val_element::Val{:iron}) = 
    group_element_tracers(nutrients.iron, bgc, val_element)
group_element_tracers(nutrient::SingleTracerNutrient, ::NPD{FT}, val_element) where FT = 
    NamedTuple{(Symbol(nutrient),)}((one(FT), ))
group_element_tracers(::NitrateAmmonia, ::NPD{FT}, val_element) where FT = 
    (NO₃ = one(FT), NH₄ = one(FT))

group_element_tracers(::Abiotic, args...) = NamedTuple()
group_element_tracers(::ImplicitProductivity, args...) = NamedTuple()
group_element_tracers(::InstantRemineralisationDetritus, args...) = NamedTuple()
group_element_tracers(::Nothing, args...) = NamedTuple()

for thing in (PhytoZoo, Detritus, DissolvedParticulate)
    for element in (:nitrogen, :iron, :phosphate)
        ratio_name = Symbol(element, :_ratio)
        @eval begin
            function group_element_tracers(group::$thing, bgc, ::Val{$(QuoteNode(element))})
                ratio = $(ratio_name)(bgc.plankton, bgc)
                names = required_biogeochemical_tracers(group)
                return NamedTuple{names}(repeat([ratio], length(names)))
            end
        end
    end
end

for thing in (PhytoZoo, Detritus)
    @eval begin
        function group_element_tracers(group::$thing, bgc, ::Val{:carbon}) # add specialisation for explicit calcite when done
            ratio = carbon_ratio(bgc.plankton, bgc)
            rain_ratio = calcite_rain_ratio(bgc.plankton, bgc)

            names = required_biogeochemical_tracers(group)

            return NamedTuple{names}(repeat([ratio * (1 + rain_ratio)], length(names)))
        end
    end
end

function group_element_tracers(::DissolvedParticulate{N, M, DN, PN}, bgc, ::Val{:carbon}) where {N, M, DN, PN} # add specialisation for explicit calcite when done
    ratio = carbon_ratio(bgc.plankton, bgc)
    rain_ratio = calcite_rain_ratio(bgc.plankton, bgc)
    
    ratios = (repeat([ratio], N)...,
              repeat([ratio * (1 + rain_ratio)], M)...)

    return NamedTuple{(DN..., PN...)}(ratios)
end

group_element_tracers(::CarbonateSystem, args...) = NamedTuple()
group_element_tracers(::CarbonateSystem, bgc::NPD{FT}, ::Val{:carbon}) where FT = 
    (; DIC = one(FT))

group_element_tracers(::Oxygen, args...) = NamedTuple()
group_element_tracers(::Oxygen, bgc::NPD{FT}, ::Val{:oxygen}) where FT = 
    (; O₂ = one(FT))

for thing in (PhytoZoo, Detritus, DissolvedParticulate)
    @eval begin
        function group_element_tracers(group::$thing, bgc::NPD{<:Any, <:Any, <:Any, <:Any, <:Any, <:Oxygen}, ::Val{:oxygen}) 
            ratio = carbon_ratio(bgc.plankton, bgc)
            rO = - bgc.oxygen.production_oxygen_carbon_ratio

            names = required_biogeochemical_tracers(group)

            return NamedTuple{names}(repeat([ratio * rO], length(names)))
        end
    end
end

group_element_tracers(::Nutrients{<:NitrateAmmonia}, bgc::NPD{<:Any, <:Any, <:Any, <:Any, <:Any, <:Oxygen}, ::Val{:oxygen}) = 
    (; NO₃ = carbon_ratio(bgc.plankton, bgc) * bgc.oxygen.nitrification_oxygen_carbon_ratio)

# got to give up on dispatch for this
function group_element_tracers(nutrients::Nutrients, bgc::NPD{<:Any, <:Any, <:Any, <:Any, <:Any, <:Oxygen}, ::Val{:oxygen})
    if nutrients.nitrogen isa SingleTracerNutrient
        return (; N = carbon_ratio(bgc.plankton, bgc) /
                      nitrogen_ratio(bgc.plankton, bgc) * bgc.oxygen.nitrification_oxygen_carbon_ratio)
    elseif nutrients.phosphate isa SingleTracerNutrient
        return (; PO₄ = carbon_ratio(bgc.plankton, bgc) /
                        phosphate_ratio(bgc.plankton, bgc) * bgc.oxygen.nitrification_oxygen_carbon_ratio)
    elseif nutrients.iron isa SingleTracerNutrient
        return (; Fe = carbon_ratio(bgc.plankton, bgc) /
                       iron_ratio(bgc.plankton, bgc) * bgc.oxygen.nitrification_oxygen_carbon_ratio)
    elseif nutrients.silicate isa SingleTracerNutrient
        return (; Si = carbon_ratio(bgc.plankton, bgc) /
                       silicon_ratio(bgc.plankton, bgc) * bgc.oxygen.nitrification_oxygen_carbon_ratio)
    else
        return NamedTuple()
    end
end

available_nutrients(nutrients) = 
    [name for name in propertynames(nutrients) if !isnothing(getproperty(nutrients, name))]

group_element_tracers(::CarbonNitrogenDissolvedParticulate{FT}, bgc, ::Val{:nitrogen}) where FT =
    (; DON = one(FT), sPON = one(FT), bPON = one(FT))

function group_element_tracers(::CarbonNitrogenDissolvedParticulate, bgc::NPD, ::Val{:phosphate})
    R_PN = phosphate_ratio(bgc.plankton, bgc) / 
          nitrogen_ratio(bgc.plankton, bgc)

    return (; DON = R_PN, sPON = R_PN, bPON = R_PN)
end

function group_element_tracers(::CarbonNitrogenDissolvedParticulate, bgc::NPD, ::Val{:iron})
    R_FeN = iron_ratio(bgc.plankton, bgc) / 
            nitrogen_ratio(bgc.plankton, bgc)

    return (; DON = R_FeN, sPON = R_FeN, bPON = R_FeN)
end

function group_element_tracers(::CarbonNitrogenDissolvedParticulate{FT}, bgc, ::Val{:carbon}) where FT
    rain_ratio = calcite_rain_ratio(bgc.plankton, bgc)

    return (; DOC = one(FT), sPOC = one(FT) + rain_ratio, bPOC = one(FT) + rain_ratio)
end

function group_element_tracers(::CarbonNitrogenDissolvedParticulate, bgc::NPD, ::Val{:oxygen})
    rO = - bgc.oxygen.production_oxygen_carbon_ratio

    return (; DOC = rO, sPOC = rO, bPOC = rO)
end

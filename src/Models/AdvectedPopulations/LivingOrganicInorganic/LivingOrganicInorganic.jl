module LivingOrganicInorganicModel

include("Living/Living.jl")
include("Organic/Organic.jl")
include("Inorganic/Inorganic.jl")

struct LivingOrganicInorganic{L, O, I}
       living :: L
      organic :: O
    inorganic :: I
end


end # module
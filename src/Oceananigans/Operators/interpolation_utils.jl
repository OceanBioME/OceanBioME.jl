using OceanBioME.Oceananigans.Fields: isacolumn

import Oceananigans.Operators: identity1, identity2, identity3, identity4, identity5, identity6, number_of_identities
#I don't really understand what the whole identity1, identity2, ... thing is about and why everything can't just use identity1? 

for i = 1:number_of_identities
    identity = Symbol(:identity, i)

    @eval begin
        @inline $identity(i, j, k, grid, c) = @inbounds isacolumn(c) ? c[i, j, 1] : c[i, j, k]
        @inline $identity(i, j, k, grid, a::Number) = a
        @inline $identity(i, j, k, grid, F::TF, args...) where TF<:Function = F(i, j, k, grid, args...)
    end
end
# this should be in Oceananigans

using Atomix, CUDA

using Oceananigans: Field
using OffsetArrays: OffsetArray

import Atomix: pointer

@inline function pointer(ref::Atomix.Internal.IndexableRef{<:Field, Tuple{Vararg{Int64, N}}} where {N})
    i = LinearIndices(ref.data.data)[ref.indices...]
    return Base.pointer(ref.data.data, i)
end

# Field on CPU
function atomic_add!(fld::Field, i, j, k, value)
    Atomix.@atomic fld[i, j, k] += value
end 

# Field on GPU which is adapted to field.data
function atomic_add!(fld::OffsetArray, i, j, k, value)
    linear_index = LinearIndices(fld)[i, j, k]
    CUDA.@atomic fld.parent[linear_index] += value
end

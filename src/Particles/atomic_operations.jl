# this should be in Oceananigans

using Atomix

using Oceananigans: Field
using OffsetArrays: OffsetArray

# Field outside of a kernel
@inline atomic_add!(fld::Field, i, j, k, value) = atomic_add!(fld.data, i, j, k, value)

# Field inside a kernel, which is adapted to field.data on every architecture.
# The atomic is applied to the parent array so that backend `Atomix` extensions
# (e.g. for `CuDeviceArray`) dispatch on the device array type
@inline function atomic_add!(fld::OffsetArray, i, j, k, value)
    linear_index = LinearIndices(fld)[i, j, k]
    Atomix.@atomic fld.parent[linear_index] += value
    return nothing
end

# this should be in Oceananigans

using Atomix

using Oceananigans: Field
using OffsetArrays: AbstractOffsetArray, offsets

# Field outside of a kernel
@inline atomic_add!(fld::Field, i, j, k, value) = atomic_add!(fld.data, i, j, k, value)

# Field inside a kernel, which is adapted to field.data on every architecture: a
# `StaticStridesOffsetArray` since Oceananigans 0.117, an `OffsetArray` before. The atomic is
# applied to the parent array so that backend `Atomix` extensions (e.g. for `CuDeviceArray`)
# dispatch on the device array type
@inline function atomic_add!(fld::AbstractOffsetArray, i, j, k, value)
    data = parent(fld)
    oi, oj, ok = offsets(fld)
    linear_index = LinearIndices(data)[i - oi, j - oj, k - ok]
    Atomix.@atomic data[linear_index] += value
    return nothing
end

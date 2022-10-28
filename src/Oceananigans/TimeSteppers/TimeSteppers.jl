module TimeSteppers

using OceanBioME.Oceananigans.Fields: AuxiliaryFields, isacolumn
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device, device_event
using Oceananigans: prognostic_fields
using Oceananigans.TurbulenceClosures: implicit_step!
using KernelAbstractions

SlicedIndices=Union{Tuple{UnitRange{Int}, Colon, Colon}, Tuple{Colon, UnitRange{Int}, Colon}, Tuple{Colon, Colon, UnitRange{Int}}, Tuple{Colon, UnitRange{Int}, UnitRange{Int}}, Tuple{UnitRange{Int}, Colon, UnitRange{Int}}, Tuple{UnitRange{Int}, UnitRange{Int}, Colon}}
Union{Tuple{Colon, Colon, UnitRange{Int64}}, Tuple{Colon, UnitRange{Int64}, Colon}, Tuple{Colon, UnitRange{Int64}, UnitRange{Int64}}, Tuple{UnitRange{Int64}, Colon, Colon}, Tuple{UnitRange{Int64}, Colon, UnitRange{Int64}}, Tuple{UnitRange{Int64}, UnitRange{Int64}, Colon}}

include("quasi_adams_bashforth_2.jl")
include("runge_kutta_3.jl")
include("store_tendencies.jl")

end

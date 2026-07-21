using OceanBioME, Test, CUDA, Oceananigans, JLD2, Oceananigans.Units, Documenter

using Oceananigans.Fields: ConstantField

if get(ENV, "CUDA_VISIBLE_DEVICES", nothing) == "-1"
    architecture = CPU()
elseif CUDA.functional()
    architecture = GPU()
else
    error("CUDA is not functional but CUDA_VISIBLE_DEVICES is not set to -1. GPU tests cannot run without a working GPU.")
end
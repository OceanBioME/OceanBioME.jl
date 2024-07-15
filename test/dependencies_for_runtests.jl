using OceanBioME, Test, CUDA, Oceananigans, JLD2#, Documenter

architecture = CUDA.has_cuda() ? GPU() : CPU()
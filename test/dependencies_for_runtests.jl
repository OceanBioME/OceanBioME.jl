using OceanBioME, Test, CUDA, Oceananigans, JLD2, Oceananigans.Units#, Documenter

architecture = CUDA.has_cuda() ? GPU() : CPU()
using Pkg

Pkg.instantiate()

Pkg.add("Test", "CUDA", "DataDeps", "Documenter", "Statistics", "JLD2")

Pkg.precompile()

using Oceananigans

arch = GPU()

include("runtests.jl")
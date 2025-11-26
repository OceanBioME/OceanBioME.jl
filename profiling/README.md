# A package doing profiling of OceanBioME

This package is hosted inside the OceanBioME repository, but in fact it just uses it as a
dependency. Hence one needs to take care since it will be easy to profile not a local version of 
OceanBioME in the parent directory, but one fetched from the registry.

To generate profile just instantiate the environment and run the script:

```julia
using Pkg
Pkg.activate("profiling")
Pkg.instantiate()
include("script.jl")
```

This should produce an interactive SVG file `prof.svg`.

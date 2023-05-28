# Visualize output

In the examples we use [Makie.jl](https://docs.makie.org/stable/) for plotting.

Makie comes with a few [backends](https://docs.makie.org/stable/#makie_ecosystem). In the documented examples
we use [CairoMakie](https://docs.makie.org/stable/documentation/backends/cairomakie/) since this backend
works well on headless devices, that is, devices without monitor. Because the documentation is automatically
built via GitHub actions the CairoMakie backend is necessary. However, users that want to run the examples on
devices with a monitor might want to change to [GLMakie](https://docs.makie.org/stable/documentation/backends/glmakie/)
that displays figures in an interactive window. To do that you need to install GLMakie, e.g.,

```julia
using Pkg
pkg"add GLMakie"
```

and replace `using CairoMakie` with `using GLMakie`.

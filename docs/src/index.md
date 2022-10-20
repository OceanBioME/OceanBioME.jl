# OceanBioME.jl

```@meta
CurrentModule = OceanBioME
DocTestSetup = quote
    using OceanBioME
end
```
```@autodocs
Modules = [OceanBioME]
```

OceanBioME.jl aims to be a flexible and friendly ocean biogeochemical modelling environment for simulating coupled physics-biogeochemical models across a range of ocean scenarios. It is built within [Oceananigans.jl](https://github.com/CliMA/Oceananigans.jl) which resolves the transport of the biogeochemical tracers and individuals and provides an excellent framework for implementing advanced features such as complex boundary conditions (e.g. sediment models).

OceanBioME.jl is currently developed primarily through the generous and greatly appreciated support of the [Center for Climate Repair at Cambridge](https://www.climaterepair.cam.ac.uk/) and the [Gordon and Betty Moore Foundation](https://www.moore.org/). 

OceanBioME.jl currently provides a core of several biogeochemical models (NPZD and [LOBSTER](https://doi.org/10.1029/2004JC002588), a medium complexity model, and [PISCES](https://doi.org/10.5194/gmd-8-2465-2015) in an early stage of testing), air-sea gas exchange models to provide appropriate top boundary conditions, and example sediment models for the benthic boundary.

As it was originally conceived for the study of submesoscale systems such as kelp forests we have also implimented a framework for integrating the growth of biological/active Lagrangian particles which move around and can interact with the (Eulerian) tracer fields - consuming nutrients and carbon dioxide, and depositing waste. Currently, included is a simple sugar kelp growth model which can be used in a variety of dynamical scenarios such as free floating or fixed to the ground, with more macroalgae models expected to be ready soon.

An overview of the current and planned future design of OceanBioME is shown in the diagram below:
![Diagram of high level structure of OceanBioME.jl showing interlink between different components ... (improve this alt text)](overview.png)

Future development ideas include:
- [Completing testing of PISCES](https://github.com/OceanBioME/OceanBioME.jl/tree/PISCES) and perform validation experiments
- Explore the potential for fixed plant populations (e.g. seagrass) models which may be implemented with the [forced column fields developed for sediment modelling](numerical_implimentation/sediments.md)
- Complete changes for [GPU compatibility](https://github.com/OceanBioME/OceanBioME.jl/pull/23) (may require upstream Oceananigans.jl or even CUDA.jl changes)
- Expand examples (e.g. global circulation example)

Please see [this project board](https://github.com/orgs/OceanBioME/projects/4) for our current roadmap/progress.


## Quick install

OceanBioME ~~is a~~ will be a [registered Julia package](https://julialang.org/packages/). So to install it,

1. [Download Julia](https://julialang.org/downloads/).

2. Launch Julia and type

```julia
julia> using Pkg

julia> Pkg.add("https://github.com/OceanBioME/OceanBioME.jl")
```
If you wish to run the examples, until we move to a better system, please then move into the examples directory and download our example data files:
```shell
shell> cd OceanBioME/examples

shell> git clone --depth 1 USERNAME@github.com:OceanBioME/OceanBioME_example_data.git
```
!!! compat "Julia 1.6 or newer"
    The latest version of OceanBioME strongly suggests _at least_ Julia 1.8 or later to run.
    While most scripts will run on Julia 1.6 or 1.7, OceanBioME is _only_ tested on Julia 1.8.

If you're [new to Julia](https://docs.julialang.org/en/v1/manual/getting-started/) and its [wonderful `Pkg` manager](https://docs.julialang.org/en/v1/stdlib/Pkg/), Oceananigans provide a very useful [wiki](https://github.com/CliMA/Oceananigans.jl/wiki) with [more detailed installation instructions](https://github.com/CliMA/Oceananigans.jl/wiki/Installation-and-getting-started-with-Oceananigans).

## Places to find OceanBioME information

* This documentation, which provides
    * example scripts,
    * explanations of model implementation methods,
    * details of currently implimented models, and
    * a library documenting all user-facing objects and functions.
* [Discussions on the OceanBioME github](https://github.com/OceanBioME/OceanBioME.jl/discussions)
  
    If you've got a question or something to talk about, don't hesitate to [start a new discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions/new?)!

* [Issues](https://github.com/OceanBioME/OceanBioME.jl//issues) and [pull requests](https://github.com/OceanBioME/OceanBioME.jl/pulls) also contain lots of information about problems we've found, solutions we're trying to implement, and ideas for the future.

## Getting in touch

Whether you need help getting started with OceanBioME, found a bug, want OceanBioME to be more expanded, or just want to chat about our project, you can:

* [Start a discussion](https://github.com/OceanBioME/OceanBioME.jl/discussions). 
* [Open an issue](https://github.com/OceanBioME/OceanBioME.jl/issues). Issues are best if you think the OceanBioME source code needs attention: a bug, a sign error, an important missing feature, or a typo in this documentation.

## Citing

If you use OceanBioME as part of your research, teaching, or other activities, we would be grateful if you could
cite our work and mention OceanBioME by name, as well as citing and [acknowledging Oceananigans](https://clima.github.io/OceananigansDocumentation/stable/#Citing) as without them this package would not be possible.

We do not currently have a citation for OceanBioME so please reach out if you wish to cite it, and we will expedite the process of [making it citable](https://joss.theoj.org/about).

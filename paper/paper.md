---
title: 'OceanBioME.jl: A flexible environment for modelling the coupled interactions between ocean biogeochemistry and physics'
tags:
  - julia
  - biogeochemistry
  - climate
  - ocean
  - carbon
authors:
  - name: Jago Strong-Wright
    orcid: 0000-0002-7174-5283
    corresponding: true
    affiliation: "1, 2"
  - name: John R Taylor
    affiliation: "1, 2"
  - name: Si Chen
    affiliation: "1, 2"
  - name: Gregory LeClaire Wagner
    orcid: 0000-0001-5317-2445
    affiliation: 3
  - name: Navid C Constantinou
    affiliation: 4, 5
  - name: Simone Silvestri
    affiliation: 3
affiliations:
 - name: Department of Applied Mathematics and Theoretical Physics, University of Cambridge, Cambridge, United Kingdom
   index: 1
 - name: Centre for Climate Repair, Cambridge, United Kingdom
   index: 2
 - name: Massachusetts Institute of Technology
   index: 3
 - name: Australian National University
   index: 4
 - name: Australian Research Council Centre for Climate Extremes
   index: 5
date: 15 March 2023
bibliography: paper.bib
---

# Summary

``OceanBioME.jl`` is a flexible modelling environment written in Julia [@julia] for modelling the coupled interactions between ocean biogeochemistry, carbonate chemistry, and physics.
``OceanBioME.jl`` can be used as a stand-alone box model, or integrated into ``Oceananigans.jl`` [@Oceananigans] simulations of ocean-flavoured fluid dynamics in one-, two-, or three-dimensions.
As a result, ``OceanBioME.jl`` and ``Oceananigans.jl`` can be used to simulate the biogeochemical response across an enormous range of scales: from surface boundary layer turbulence at the meter scale to eddying global ocean simulations at the planetary scale, and on computational systems ranging from laptops to supercomputers.
An example of a problem involving small-scale flow features is shown in \autoref{fig1}.
``OceanBioME.jl`` leverages Julia's multiple dispatch and effective inline capabilities to fuse its computations directly into existing ``Oceananigans.jl`` kernels, thus maintaining ``Oceananigans.jl``'s bespoke performance, memory- and cost-efficiency on GPUs in ``OceanBioME.jl``-augmented simulations.

![Here we replicate the Eady problem where a background buoyancy gradient and corresponding thermal wind generate a sub-mesoscale eddy, roughly following the setup of Taylor (2016).
To this physical setup, we added a medium complexity (9 tracers) biogeochemical model, some of which are shown above.
On top of this, we added particles modelling the growth of sugar kelp which are free-floating and advected by the flow, and carbon dioxide exchange from the air.
Thanks to Julia's speed and efficiency the above model (1 km × 1 km × 100 m with 64 × 64 × 16 grid points) took about 30 minutes of computing time to simulate 10 days of evolution on an Nvidia P100 GPU. Figure made with `Makie.jl` [@makie]. \label{fig1}](eady_example.png)

``OceanBioME.jl`` is built with a highly modular design that allows user control and customization.
There are two distinct module types implemented in ``OceanBioME.jl``:

- First, tracer-based ecosystem modules are formulated in `AdvectedPopulations` as a set of coupled ordinary differential equations.
These equations can be solved by ``OceanBioME.jl`` as box models, which is particularly useful for testing.
The same modules can be integrated by ``Oceananigans.jl`` to provide tracer-based ecosystem models.

- The second module type is `Individual` "biologically active" particles.
These consist of individual-based models which are solved along particle paths and can be coupled with the tracer-based modules and physics from ``Oceananigans.jl``.
The biologically active particles can be advected by the currents, and/or they can move according to prescribed dynamics.
For example, migrating zooplankton or fish can be modelled with biologically active particles and ``OceanBioME.jl`` allows these to interact with tracer-based components such as phytoplankton or detritus.

`AdvectedPopulations` are supported by `Boundaries` modules which provide information at the top and bottom of the ocean. For example, the GasExchange submodule calculate the flux of carbon dioxide and oxygen at the sea surface, while sediment modules calculate fluxes of carbon and oxygen at the seafloor.

We provide a simple framework and utilities (such as light attenuation integration) to build the necessary components of biogeochemical models.
With the provided models, currently a simple Nutrient-Phytoplankton-Zooplankton-Detritus [@npzd] model, and an intermediate complexity model, LOBSTER [@lobster], we have set up a straightforward "plug and play" framework to add additional tracers such as carbonate and oxygen chemistry systems, and additional forcing.
Additionally, we have implemented comprehensive air-sea flux models [e.g. @wanninkhof:1992] and sediment models [e.g. @soetaert:2000] which can easily be applied to tracers in the models.
We focus on the simulation of idealized sub-mesoscale systems, but this flexible framework allows users to model problems of any scale.
This framework is made possible by our contributions to ``Oceananigans.jl``, adding a streamlined user interface to swap biogeochemical models with no modification to other model configurations.
This interface also facilitates rapid prototyping, as models can be implemented and swapped easily by just extending a few key functions.
This flexibility and ease-of-use is unmatched in existing biogeochemical models.

``OceanBioME.jl`` was designed specifically to study ocean carbon dioxide removal (OCDR) strategies.
Assessing the effectiveness and impacts of OCDR is challenging due to the complexities of the interactions between the biological, chemical, and physical processes involved in the carbon cycle.
Moreover, field trials of OCDR interventions are generally small-scale and targeted, while the intervention required to have a climate-scale impact is regional or global.
We have built ``OceanBioME.jl`` to meet these challenges by creating tools that provide a modular interface to the different components within the ocean modelling framework provided by ``Oceananigans.jl``.
This allows easy access to a suite of biogeochemical models ranging from simple idealized to full-complexity models.
\autoref{fig2} shows a simple column model with an OCDR intervention (macroalgae growth) added after a warm-up period, which increases the carbon export of the system.

![Here we show the results of a 1D model, forced by idealised light and mixing, which qualitatively reproduces the biogeochemical cycles in the North Atlantic.
We then add kelp (500 frond / m² in the top 50 m of water) in December of the 2ⁿᵈ year (black vertical line) which causes an increase in air-sea carbon dioxide exchange and sinking export, as well as a change in the phytoplankton growth cycle.
Figure made with `Makie.jl` [@makie]. \label{fig2}](column_example.png)

The biologically active particles built into ``OceanBioME.jl`` are particularly useful for OCDR applications.
Accurate carbon accounting is essential for assessing the effectiveness of OCDR strategies.
Biologically active particles can be used to track carbon from a particular source while accounting for interactions with its surroundings.
Biologically active particles can also be used to model OCDR deployment strategies including seaweed cultivation, alkalinity enhancement, and marine biomass regeneration.
``OceanBioME.jl`` currently includes an extended version of the sugar kelp model presented by @broch:2012 as an example of the utility and implementation of these features.

The implementation of OceanBioME.jl models allows for seamless integration with data assimilation packages, such as ``EnsembleKalmanProcesses.jl`` [@ekp]. This feature facilitates rapid calibration of model parameters, providing a powerful utility for integrating observations and models, with the potential to improve model skill and identify key sources of uncertainty.

A key metric for the validity of biogeochemical systems is the conservation of elements such as carbon and nitrogen in the system.
We therefore continuously test the implemented models in a variety of simple scenarios (i.e. isolated, with/without air-sea flux, with/without sediment) to ensure basic conservations are fulfilled, and will continue to add tests for any new models.
Additionally, we check ``OceanBioME.jl`` utilities through standard tests such as comparison to analytical solutions for light attenuation, and conservation of tracers for active particle exudation and sinking.

Finally, this software is currently facilitating multiple research projects into ocean CDR which would have been significantly harder with other solutions.
For example, Chen (In prep.) is using the active particle coupling provided to investigate the effects of location and planting density of kelp in the open ocean on their carbon drawdown effect, as in the example above.
Additionally, Strong-Wright (In prep.) is using the coupling of both the biogeochemistry and easy interface to couple the physics to study flow interactions with a fully resolved giant kelp forest model including the effects on nutrient transport and distribution.

# Acknowledgements

We would like to thank the [Climate Modeling Alliance](https://clima.caltech.edu) team and ``Oceananigans.jl`` contributors for their fantastic project. We are also very grateful for the support and funding of the [Centre for Climate Repair at Cambridge](https://www.climaterepair.cam.ac.uk/) and the [Gordon and Betty Moore Foundation](https://www.moore.org/).

# References

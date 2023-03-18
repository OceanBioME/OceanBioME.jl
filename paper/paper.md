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
affiliations:
 - name: Department of Applied Mathematics and Theoretical Physics, University of Cambridge, Cambridge, United Kingdom
   index: 1
 - name: Centre for Climate Repair at Cambridge, Cambridge, United Kingdom
   index: 2
date: 15 March 2023
bibliography: paper.bib
---

# Summary

Modelling the oceans biological systems is challenging due to the complex interplay of physical, chemical, geological, and biological factors, the mechanisms of which may be poorly understood. An area where improvements in this area are vital is the study of the effectiveness and impacts of ocean based carbon dioxide removal (CDR) strategies. With this application in mind we have written ``OceanBioME.jl``, an ocean biogeochemical modelling environment intended to allow easy access to flexible ocean biogeochemistry modelling with the easy addition of new subsystems. 

In ``OceanBioME.jl`` in collaboration with ``Oceananigans.jl`` [@Oceananigans] we provide a simple framework and utilities (such as light attenuation integration) to couple all the necessary components of biogeochemical models. With the provided models, currently a simple NPZD [@npzd] model, an intermediate complexity model LOBSTER [@lobster], and PISCES [@pisces] we have set up a straightforward "plug and play" framework to add additional tracers such as carbonate and oxygen chemistry systems, and additional forcing. Additionally, we have implemented a comprehensive air-sea flux models [e.g. @wanninkhof:1992] and experimentally sediment models [e.g. @soetaert:2000] which can easily be applied to tracers in the models. We focus on the simulation of idealized sub-mesoscale systems, but this flexible framework allows users to model problems of any scale.

In order to simulate the effects of additional systems such as macroalgae we have expanded upon ``Oceananigans.jl`` *Lagrangian particles* to provide a framework for "active" particles where their internal dynamics can interact with tracer fields (e.g. nutrient uptake). We also include an extended version of the sugar kelp model presented in @broch:2012 as an example of the utility and implimentation of these features. 

We have also formulated the models such that they are easy to use alongside data assimilation packages such as ``EnsembleKalmanProcesses.jl`` [@ekp] to calibrate their parameters. This provides vital utility for integrating observations and will allow improved validation of CDR strategies.

![Here we show the results of a 1D column model, forced by idealised light and mixing, forced with a seasonal cycle of temperature and light, which qualitatively reproduces the biogeochemical cycles in the north Atlantic. We then add kelp in December of the 2nd year which causes an increase in air-sea carbon dioxide exchange and sinking export, as well as a change in the phytoplankton growth cycle. Plots made with `Makie` [@makie].](column_example.png)

![Here we replicate the Eady problem where a background buoyancy gradient and corresponding thermal wind generate submesoscale eddies, roughly following the setup of Taylor (2016). To this physical setup we added a medium complexity (9 tracers) biogeochemical model, some of which are shown above. On top of this we added particles modelling the growth of sugar kelp which are free-floating and advected by the flow, and carbon dioxide exchange from the air. A key advantage of writing this package in Julia is it offers accessibility similar to highl-level languages such as Python, with the speed of languages like C and Fortran and built-in parallelism. This means that models can be run significantly faster than the equivalent in other high-level languages. Additionally OceanBioME can run on GPUs, allowing the above model (1km x 1km x 100m with 64 x 64 x 16 grid points) to simulate 10 days of evolution in about 30 minutes of computing time.](eady_example.png)

A key metric for the validity of biogeochemical systems is the conservation of nitrogen in the system. We therefore continuously test the implemented models in a variety of simple scenarios (i.e. isolated, with/without air-sea flux, with/without sediment) to ensure basic conservations are fulfilled, and will continue to add tests for any new models. Additionally, we ensure the validity of the utilities provided through unit tests such as comparison to analytical solutions for light attenuation, and conservation of tracers for active particle exudation and sinking.

Flexible biogeochemical modelling frameworks similar to ``OceanBioME.jl`` are uncommon and tend to require more significant knowledge of each coupled system, a more cumbersome configuration process, provide a narrower breadth of utility, are not openly available, or are more computationally intensive. For example among the open source alternatives NEMO [@nemo] provides a comprehensive global biogeochemical modelling framework but requires complex configuration and is unsuited for local ecosystem modelling, while MACMODS [@macmods] provides more limited functionality on a slower platform.

Finally, this software is currently facilitating multiple research projects into ocean CDR which would have been significantly harder with other solutions. For example, Chen (In prep.) is using the active particles coupling provided to investigate the effects of location and planting density of kelp in the open ocean on their carbon drawdown effect, as in the example above. Additionally, Strong-Wright (In prep.) is using the coupling of both the biogeochemistry and easy interface to couple the physics to study flow interactions with a fully resolved giant kelp forest including the effects on nutrient transport and distribution.

[comment]: <> (Not convinved we need this section since Oceananigans doesn't have one, the above is already about the same length as their paper, and it doesn't really flow)
# Statement of need

The UN declared 2021-2030 the *Decade of Ocean Science for Sustainable Development*, among the goals identified by this declaration are "the protection and restoration of ecosystems and biodiversity" and "unlocking ocean-based solutions to climate change" while being "transparent and accessible" to empower all stakeholders in their decision-making [@un:2018]. Improving our quantitative understanding of ocean ecosystems must be a priority in achieving these goals to allow their informed management and adaptation [@ryabinin:2019], and ``OceanBioME.jl`` contributes to these goals by facilitating easier and more comprehensive access to ecosystem modelling. 

# Acknowledgements

We would like to thank the ``Oceananigans`` contributors for their fantastic project, and particularly Greg Wanger for his advice and support, we are also very grateful for the support and funding of the [Centre for Climate Repair at Cambridge](https://www.climaterepair.cam.ac.uk/) and the [Gordon and Betty Moore Foundation](https://www.moore.org/).

# References

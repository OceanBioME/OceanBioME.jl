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
    equal-contrib: true
    corresponding: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: John R Taylor
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1, 2"
  - name: Si Chen
    affiliation: 1
affiliations:
 - name: Department of Applied Mathematics and Theoretical Physics, University of Cambridge, Cambridge, United Kingdom
   index: 1
 - name: Centre for Climate Repair at Cambridge, Cambridge, United Kingdom
   index: 2
date: 18 August 2020
bibliography: paper.bib
---

# Summary

Modelling the oceans biological systems is challenging due to the complex interplay of physical, chemical, geological, and biological factors, the mechanisms of which may be poorly understood. An area where improvements in this area are vital is the study of the effectiveness and impacts of ocean based carbon dioxide removal (CDR) strategies. With this application in mind we have written ``OceanBioME.jl``, an ocean biogeochemical modelling environment intended to allow easy access to flexible complexity ocean biogeochemistry modelling with the easy addition of new subsystems. 

In ``OceanBioME.jl`` we provide a framework and utilities (such as light attenuation integration) to couple biogeochemical models [e.g. LOBSTER @lobster] with the physics model ``Oceananigans.jl`` [@oceananigans], implemented as forced tracer fields. With the provided models (currently LOBSTER and a simple NPZ model) we have set up the framework to be straightforward to turn on and off additional tracers such as carbonate and oxygen chemistry systems. Additionally, we have implemented various air-sea flux models [e.g. @wanninkhof:1992] and sediment models [e.g. @soetaert:2000] which can easily be applied to arbitrary tracers in the models. 

In order to simulate the effects of additional systems such as macroalgae we have expanded upon ``Oceananigans.jl`` *Lagrangian particles* to provide a framework for "active" particles where their internal dynamics and initial positions are provided by the user, along with their field dependencies and "source/sink" variables (e.g. nutrient uptake). The particles are then integrated, and interact, with the tracers. We include an extended version of the sugar kelp model presented in @broch:2012 as an example of the utility of this feature.

![Fig. 1](example.png)
Fig. 1: (Left) Replication of phytoplankton seasonal cycles in the subpolar regions showing the responses to changing in mixing layer depth, and light and nutrient availability, including deep spring blooms and later summer deep chlorophyll maxima. Additionally, the model is configured with a large amount of sugar kelp added in the third year showing a deepening of the phytoplankton response due to nutrient redistribution and changes to light attenuation. (Right) The air-sea flux of carbon dioxide provides useful insight into the effects of changes to the ecosystem on its CDR, here showing an increase in downward carbon flux after the kelp is added. Plots made with `GLMakie` [@glmakie].

A key metric for the validity of biogeochemical systems is the conservation of nitrogen in the system. We therefore continuously test the implemented models in a variety of simple scenarios (i.e. isolated, with/without air-sea flux, with/without sediment) to ensure basic conservations are fulfilled, and will continue to add tests for any new models. Additionally, we ensure the validity of the utilities provided through unit tests such as comparison to analytical solutions for light attenuation, and conservation of tracers for active particle exudation and sinking.

Flexible biogeochemical modelling frameworks similar to ``OceanBioME.jl`` are uncommon and tend to require more significant knowledge of each coupled system, a more cumbersome configuration process, provide a narrower breadth of utility, are not openly available, or are more computationally intensive. For example among the open source alternatives NEMO [@nemo] provides a comprehensive global biogeochemical modelling framework but requires complex configuration and is unsuited for local ecosystem modelling, while MACMODS [@macmods] provides more limited functionality on a slower platform.

[comment]: <> (Not convinved we need this section since Oceananigans doesn't have one, the above is already about the same length as their paper, and it doesn't really flow)
# Statement of need

The UN declared 2021-2030 the *Decade of Ocean Science for Sustainable Development*, among the goals identified by this declaration are "the protection and restoration of ecosystems and biodiversity" and "unlocking ocean-based solutions to climate change" while being "transparent and accessible" to empower all stakeholders in their decision-making [@un:2020]. Improving our quantitative understanding of ocean ecosystems must be a priority in achieving these goals to allow their informed management and adaptation [@ryabinin:2019], and ``OceanBioME.jl`` contributes to these goals by facilitating easier and more comprehensive access to ecosystem modelling. 

# Acknowledgements

We are very grateful for the support and funding of the [Centre for Climate Repair at Cambridge](https://www.climaterepair.cam.ac.uk/) and the [Gordon and Betty Moore Foundation](https://www.moore.org/). And would like to thank the ``Oceananigans.jl`` team for their fantastic project and continuous advice throughout.

# References
# [PISCES (Pelagic Interactions Scheme for Carbon and Ecosystem Studies) model](@id PISCES)
PISCES ([PEES-kays, /ˈpiːs.keːs/](https://forvo.com/word/pisc%C4%93s/#la)) is a high complexity ocean biogeochemical model with 24 prognostic tracers. 
It has previously been used with the [NEMO](https://www.nemo-ocean.eu/) transport model in the [IPSL-CM5A-LR](https://doi.org/10.1007/s00382-012-1636-1) and [CNRM-CM5](https://doi.org/10.1007/s00382-011-1259-y) CMIP-5 earth system models (ESM).
This is an early attempt to implement PISCES for use as a test bed in a more flexible environment to allow rapid prototyping and testing of new parametrisations as well as use in idealised experiments, additionally we want to be able to replicate the dynamics of the operational model for possible future use in a Julia based ESM as the ecosystem matures.

An overview of the model structure is available from the [PISCES community website](https://www.pisces-community.org):

![PISCES model structure](https://www.pisces-community.org/wp-content/uploads/2021/12/PISCES_Operational-1.png)

More documentation to follow..., for now see our list of key differences from [Aumont2015](@citet) and queries we have about the parametrisations on the [notable differences](@ref PISCES_queries).

## Model conservation
When the permanent scavenging of iron, and nitrogen fixation are turned off, PISCES conserves:

- Carbon: ``\partial_tP + \partial_tD + \partial_tZ + \partial_tM + \partial_tDOC + \partial_tPOC + \partial_tGOC + \partial_tDIC + \partial_tCaCO_3=0``
- Iron: ``\partial_tPFe + \partial_tDFe + \theta^{Fe}\left(\partial_tZ + \partial_tM + \partial_tDOC\right) + \partial_tSFe + \partial_tBFe + \partial_tFe=0``
- Phosphate: ``\theta^P\left(\partial_tPFe + \partial_tDFe + \partial_tZ + \partial_tM + \partial_tDOC + \partial_tPOC + \partial_tGOC\right) + \partial_tPO_4=0``
- Silicon: ``\partial_tDSi + \partial_tPSi + \partial_tSi=0``
- Nitrogen: ``\theta^N\left(\partial_tPFe + \partial_tDFe + \partial_tZ + \partial_tM + \partial_tDOC + \partial_tPOC + \partial_tGOC\right) + \partial_tNH_4 + \partial_tNO_3=0``
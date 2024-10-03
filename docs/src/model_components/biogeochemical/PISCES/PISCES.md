# [PISCES (Pelagic Interactions Scheme for Carbon and Ecosystem Studies) model](@id PISCES)
PISCES ([PEES-kays, /ˈpiːs.keːs/](https://forvo.com/word/pisc%C4%93s/#la)) is a high complexity ocean biogeochemical model with 24 prognostic tracers. 
It has previously been used with the [NEMO](https://www.nemo-ocean.eu/) transport model in the [IPSL-CM5A-LR](https://doi.org/10.1007/s00382-012-1636-1) and [CNRM-CM5](https://doi.org/10.1007/s00382-011-1259-y) CMIP-5 earth system models (ESM).
This is an early attempt to implement PISCES for use as a test bed in a more flexible environment to allow rapid prototyping and testing of new parametrisations as well as use in idealised experiments, additionally we want to be able to replicate the dynamics of the operational model for possible future use in a Julia based ESM as the ecosystem matures.

An overview of the model structure is available from the [PISCES community website](https://www.pisces-community.org):

![PISCES model structure](https://www.pisces-community.org/wp-content/uploads/2021/12/PISCES_Operational-1.png)

The default configuration of PISCES in OceanBioME is the operational/standard version with 24 tracers and can be set up by writing:

```jldoctest; filter = r".*@ OceanBioME.Models.PISCESModel*"
julia> using OceanBioME, Oceananigans

julia> grid = RectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1));

julia> biogeochemistry = PISCES(; grid)
┌ Warning: This implementation of PISCES is in early development and has not yet been validated against the operational version
└ @ OceanBioME.Models.PISCESModel ~/Documents/Projects/OceanBioME.jl/src/Models/AdvectedPopulations/PISCES/PISCES.jl:344
PISCES biogeochemical model (24 tracers) 
 Light attenuation: Multi band light attenuation model with 3 bands (:PAR₁, :PAR₂, :PAR₃)
 Sediment: Nothing
 Particles: Nothing
 Modifiers: Nothing

julia> Oceananigans.Biogeochemistry.required_biogeochemical_tracers(biogeochemistry)
(:P, :PChl, :PFe, :D, :DChl, :DFe, :DSi, :Z, :M, :DOC, :POC, :GOC, :SFe, :BFe, :PSi, :CaCO₃, :NO₃, :NH₄, :PO₄, :Fe, :Si, :DIC, :Alk, :O₂, :T, :S)

```

The parametrisations can easily be switched out when the `biogeochemistry` is constructed by setting the key word parameter, see the [API documentation](@ref library_api), although we currently do not have any of the other configurations implemented. Note that `PISCES-simple` is very similar to [`LOBSTER`](@ref LOBSTER) if that is what you are looking for.

More documentation will follow but for now the equations can be found in [Aumont2015](@citet) read along side our notes [here](@ref PISCES_queries).

## Model conservation
When the permanent scavenging of iron, nitrogen fixation, and particle sinking are turned off, PISCES conserves:

- Carbon: ``\partial_tP + \partial_tD + \partial_tZ + \partial_tM + \partial_tDOC + \partial_tPOC + \partial_tGOC + \partial_tDIC + \partial_tCaCO_3=0``
- Iron: ``\partial_tPFe + \partial_tDFe + \theta^{Fe}\left(\partial_tZ + \partial_tM + \partial_tDOC\right) + \partial_tSFe + \partial_tBFe + \partial_tFe=0``
- Phosphate: ``\theta^P\left(\partial_tPFe + \partial_tDFe + \partial_tZ + \partial_tM + \partial_tDOC + \partial_tPOC + \partial_tGOC\right) + \partial_tPO_4=0``
- Silicon: ``\partial_tDSi + \partial_tPSi + \partial_tSi=0``
- Nitrogen: ``\theta^N\left(\partial_tPFe + \partial_tDFe + \partial_tZ + \partial_tM + \partial_tDOC + \partial_tPOC + \partial_tGOC\right) + \partial_tNH_4 + \partial_tNO_3=0``
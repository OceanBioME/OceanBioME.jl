"""
    CarbonateSystem(replicates = 1)

The inorganic carbon component for the `inorganic_carbon` slot of a
[`NutrientsPlanktonDetritus`](@ref) model. It adds dissolved inorganic carbon (`DIC`) and alkalinity
(`Alk`) tracers whose evolution is driven by primary production, remineralisation of organic waste,
and calcium carbonate production/dissolution (implicit calcium carbonate: calcium carbonate is not tracked as its own tracer).

`replicates` adds independent copies of the carbonate tracers, which is useful for ensemble or
perturbation experiments; each replicate evolves with the same tendency from its own state. The
replicates are either numbered or labelled:

```julia
CarbonateSystem()                      # DIC, Alk
CarbonateSystem(2)                     # DIC1, DIC2, Alk1, Alk2
CarbonateSystem(:counterfactual)       # DIC_counterfactual, Alk_counterfactual
CarbonateSystem((:control, :release))  # DIC_control, DIC_release, Alk_control, Alk_release
```

[`carbonate_replicates`](@ref) returns the tracer names of each replicate.
"""
struct CarbonateSystem{L} <: AbstractInorganicCarbon{L, (:DIC, :Alk)} end

CarbonateSystem(replicates = 1) = CarbonateSystem{replicate_labels(replicates)}()

required_biogeochemical_auxiliary_fields(::CarbonateSystem) = tuple()

function Base.show(io::IO, c::CarbonateSystem{L}) where L
    msg = "CarbonateSystem $(required_biogeochemical_tracers(c))"

    if length(L) > 1
         msg *= "\n└── Realisations: $(length(L))"
    end

    print(io, msg)

    return nothing
end

using Oceananigans: fields, Simulation
using KernelAbstractions
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device, architecture, on_architecture
using Oceananigans.Biogeochemistry: AbstractBiogeochemistry

import Adapt: adapt_structure, adapt
import Base: summary, show
import Oceananigans.Biogeochemistry: update_tendencies!, update_biogeochemical_state!

"""
    ZeroNegativeTracers(; exclude = ())

Construct a modifier that zeroes any negative tracers excluding those listed in `exclude`.

!!! danger "Tracer conservation"
    This method is _not_ recommended as a way to preserve positivity of tracers since
    it does not conserve the total tracer.
"""
@kwdef struct ZeroNegativeTracers{E}
    exclude :: E = ()
end

function update_biogeochemical_state!(model, zero::ZeroNegativeTracers)
    for (tracer_name, tracer) in pairs(model.tracers)
        if !(tracer_name in zero.exclude)
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end

#####
##### Infastructure to rescale negative values
#####

struct ScaleNegativeTracers{FA, SA, FV, W}
           tracers :: FA
      scalefactors :: SA
invalid_fill_value :: FV
              warn :: W

    ScaleNegativeTracers(tracers::FA, scalefactors::SA, invalid_fill_value::FV, warn::W) where {FA, SA, W, FV} =
        warn ? error("Warning not currently implemented") : new{FA, SA, FV, W}(tracers, scalefactors, invalid_fill_value, warn)
end

adapt_structure(to, snt::ScaleNegativeTracers) = ScaleNegativeTracers(adapt(to, snt.tracers),
                                                                      adapt(to, snt.scalefactors),
                                                                      adapt(to, snt.invalid_fill_value),
                                                                      adapt(to, snt.warn))

"""
    ScaleNegativeTracers(; tracers, scalefactors = ones(length(tracers)), warn = false, invalid_fill_value = NaN)

Constructs a modifier to scale `tracers` so that none are negative. Use like:
```julia
modifier = ScaleNegativeTracers((:P, :Z, :N))
biogeochemistry = Biogeochemistry(...; modifier)
```
This method is better, though still imperfect, method to prevent numerical errors that lead to
negative tracer values compared to [`ZeroNegativeTracers`](@ref). Please see [discussion in
github](https://github.com/OceanBioME/OceanBioME.jl/discussions/48).

Future plans include implement a positivity-preserving timestepping scheme as the ideal alternative.

~~If `warn` is true then scaling will raise a warning.~~

`invalid_fill_value` specifies the value to set the total cell content to if the total is less than 0
(meaning that total tracer conservation can not be enforced). If the value is set to anything other than
`NaN` this scheme no longer conserves mass. While this may be useful to prevent spurious numerics leading
to crashing care should be taken that the mass doesn't deviate too much.

This scheme is similar to that used by [NEMO-PISCES](https://www.nemo-ocean.eu/), although they scale the 
tendency rather than the value, while other Earth system models simply set negative tracers to zero, for 
example [NCAR's MARBL](https://marbl-ecosys.github.io/versions/latest_release/index.html) and
[NEMO-TOPAZ2](https://zenodo.org/records/2648099), which does not conserve mass. More complicated schemes 
exist, for example [ROMS-BECS](https://zenodo.org/records/3988618) uses an implicite-itterative 
approach where each component is updated in sequence to garantee mass conservation, possibly at the 
expense of numerical precision.
"""
function ScaleNegativeTracers(tracers; scalefactors = ones(length(tracers)), invalid_fill_value = NaN, warn = false)
    if length(scalefactors) != length(tracers)
        error("Incorrect number of scale factors provided")
    end

    return ScaleNegativeTracers(tracers, scalefactors, invalid_fill_value, warn)
end

"""
    ScaleNegativeTracers(model::UnderlyingBiogeochemicalModel; warn = false)

Construct a modifier to scale the conserved tracers in `model`.

If `warn` is true then scaling will raise a warning.
"""
function ScaleNegativeTracers(bgc::AbstractBiogeochemistry, grid; invalid_fill_value = NaN, warn = false)
    tracers = conserved_tracers(bgc)
    scalefactors = on_architecture(architecture(grid), ones(length(tracers)))

    return ScaleNegativeTracers(tracers, scalefactors, invalid_fill_value, warn)
end

summary(scaler::ScaleNegativeTracers) = string("Mass conserving negative scaling of $(scaler.tracers)")
show(io::IO, scaler::ScaleNegativeTracers) = print(io, string(summary(scaler), "\n",
                                                          "└── Scalefactors: $(scaler.scalefactors)"))

function update_biogeochemical_state!(model, scale::ScaleNegativeTracers)
    workgroup, worksize = work_layout(model.grid, :xyz)

    dev = device(model.grid.architecture)

    scale_for_negs_kernel! = scale_for_negs!(dev, workgroup, worksize)

    tracers_to_scale = Tuple(model.tracers[tracer_name] for tracer_name in keys(scale.tracers))

    scale_for_negs_kernel!(tracers_to_scale, scale.scalefactors, scale.invalid_fill_value)
end

@kernel function scale_for_negs!(tracers, scalefactors, invalid_fill_value)
    i, j, k = @index(Global, NTuple)

    t, p = 0.0, 0.0

    for (idx, tracer) in enumerate(tracers)
        value = @inbounds tracer[i, j, k]
        scalefactor = @inbounds scalefactors[idx]

        t += value * scalefactor
        if value > 0
            p += value * scalefactor
        end
    end 

    t < 0 && (t = invalid_fill_value)

    for tracer in tracers
        value = @inbounds tracer[i, j, k]
        
        if value > 0
            value *= t / p
        else
            value = 0
        end

        @inbounds tracer[i, j, k] = value
    end
end
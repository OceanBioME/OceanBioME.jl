using Oceananigans: fields, Simulation
using KernelAbstractions
using KernelAbstractions.Extras.LoopInfo: @unroll
using Oceananigans.Utils: work_layout
using Oceananigans.Architectures: device, architecture, arch_array
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
    @unroll for (tracer_name, tracer) in pairs(model.tracers)
        if !(tracer_name in zero.exclude)
            parent(tracer) .= max.(0.0, parent(tracer))
        end
    end
end

#####
##### Infastructure to rescale negative values
#####

struct ScaleNegativeTracers{FA, SA, W}
    tracers :: FA
    scalefactors :: SA
    warn :: W

    ScaleNegativeTracers(tracers::FA, scalefactors::SA, warn::W) where {FA, SA, W} =
        warn ? error("Warning not currently implemented") : new{FA, SA, W}(tracers, scalefactors, warn)
end

adapt_structure(to, snt::ScaleNegativeTracers) = ScaleNegativeTracers(adapt(to, snt.tracers),
                                                                      adapt(to, snt.scalefactors),
                                                                      adapt(to, snt.warn))

"""
    ScaleNegativeTracers(; tracers, scalefactors = ones(length(tracers)), warn = false)

Constructs a modifier to scale `tracers` so that none are negative. Use like:
```julia
modifier = ScaleNegativeTracers((:P, :Z, :N))
biogeochemistry = Biogeochemistry(...; modifier)
```
This method is better, though still imperfect, method to prevent numerical errors that lead to
negative tracer values compared to [`ZeroNegativeTracers`](@ref). Please see [discussion in
github](https://github.com/OceanBioME/OceanBioME.jl/discussions/48).

Future plans include implement a positivity-preserving timestepping scheme as the ideal alternative.

If `warn` is true then scaling will raise a warning.
"""
function ScaleNegativeTracers(tracers; scalefactors = ones(length(tracers)), warn = false)
    if length(scalefactors) != length(tracers)
        error("Incorrect number of scale factors provided")
    end

    return ScaleNegativeTracers(tracers, scalefactors, warn)
end

"""
    ScaleNegativeTracers(model::UnderlyingBiogeochemicalModel; warn = false)

Construct a modifier to scale the conserved tracers in `model`.

If `warn` is true then scaling will raise a warning.
"""
function ScaleNegativeTracers(bgc::AbstractBiogeochemistry, grid; warn = false)
    tracers = conserved_tracers(bgc)
    scalefactors = arch_array(architecture(grid), ones(length(tracers)))

    return ScaleNegativeTracers(tracers, scalefactors, warn)
end

summary(scaler::ScaleNegativeTracers) = string("Mass conserving negative scaling of $(scaler.tracers)")
show(io::IO, scaler::ScaleNegativeTracers) = print(io, string(summary(scaler), "\n",
                                                          "└── Scalefactors: $(scaler.scalefactors)"))

function update_biogeochemical_state!(model, scale::ScaleNegativeTracers)
    workgroup, worksize = work_layout(model.grid, :xyz)

    dev = device(model.grid.architecture)

    scale_for_negs_kernel! = scale_for_negs!(dev, workgroup, worksize)

    tracers_to_scale = Tuple(model.tracers[tracer_name] for tracer_name in keys(scale.tracers))

    scale_for_negs_kernel!(tracers_to_scale, scale.scalefactors)
end

@kernel function scale_for_negs!(tracers, scalefactors)
    i, j, k = @index(Global, NTuple)

    t, p = 0.0, 0.0

    @unroll for (idx, tracer) in enumerate(tracers)
        value = @inbounds tracer[i, j, k]
        scalefactor = @inbounds scalefactors[idx]

        t += value * scalefactor
        if value > 0
            p += value * scalefactor
        end
    end 

    t < 0 && (t = NaN)

    @unroll for tracer in tracers
        value = @inbounds tracer[i, j, k]
        
        if value > 0
            value *= t / p
        else
            value = 0
        end

        @inbounds tracer[i, j, k] = value
    end
end
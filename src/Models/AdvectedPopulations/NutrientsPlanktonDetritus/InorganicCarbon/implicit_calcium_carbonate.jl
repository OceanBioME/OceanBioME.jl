"""
    CarbonateSystem(replicates = 1)

The inorganic carbon component for the `inorganic_carbon` slot of a
[`NutrientsPlanktonDetritus`](@ref) model. It adds dissolved inorganic carbon (`DIC`) and alkalinity
(`Alk`) tracers whose evolution is driven by primary production, remineralisation of organic waste,
and calcium carbonate production/dissolution (implicit calcium carbonate: calcium carbonate is not tracked as its own tracer).

Passing `replicates > 1` adds `replicates` independent copies of the carbonate tracers
(`DIC1`, `Alk1`, `DIC2`, …), which is useful for ensemble or perturbation experiments; each replicate
evolves with the same tendency as the base `DIC`/`Alk`.
"""
struct CarbonateSystem{N} <: AbstractInorganicCarbon{N} end

CarbonateSystem(replicates = 1) = CarbonateSystem{replicates}()

required_biogeochemical_tracers(::CarbonateSystem{1}) = (:DIC, :Alk)
required_biogeochemical_tracers(::CarbonateSystem{N}) where N = (map(n->Symbol(:DIC, n), 1:N)..., map(n->Symbol(:Alk, n), 1:N)...)
required_biogeochemical_auxiliary_fields(::CarbonateSystem) = tuple()

const NPD_DIC_Alk{FT} = NPD{FT, <:Any, <:Any, <:Any, <:CarbonateSystem}

# With `replicates > 1` the tracers are `DIC1`, `Alk1`, `DIC2`, … rather than `DIC` and `Alk`, and
# each replicate evolves with the tendency of the base tracer. Which replicate names exist is only
# fixed when the component is constructed, so they are matched against `N` here at compile time
# rather than having one method per name (see `component_tendency` in nutrients_plankton_detritus.jl).
@inline @generated function component_tendency(i, j, k, grid, ::CarbonateSystem{N}, ::Val{name},
                                               bgc::NPD_DIC_Alk{FT}, clock, fields, auxiliary_fields) where {N, name, FT}
    N > 1 || return :(zero($FT))

    for n in 1:N
        name == Symbol(:DIC, n) && return :($(Expr(:meta, :inline)); bgc(i, j, k, grid, Val(:DIC), clock, fields, auxiliary_fields))
        name == Symbol(:Alk, n) && return :($(Expr(:meta, :inline)); bgc(i, j, k, grid, Val(:Alk), clock, fields, auxiliary_fields))
    end

    return :(zero($FT))
end

Base.summary(carbonates::CarbonateSystem{1}) = 
    string("CarbonateSystem $(required_biogeochemical_tracers(carbonates))")

Base.summary(carbonates::CarbonateSystem{N}) where N = 
    string("CarbonateSystem{realisations = $N} $(required_biogeochemical_tracers(carbonates))")

function Base.show(io::IO, c::CarbonateSystem{N}) where N
    msg = "CarbonateSystem $(required_biogeochemical_tracers(c))"

    if N>1
         msg *= "\n└── Realisations: $N"
    end

    print(io, msg)

    return nothing
end

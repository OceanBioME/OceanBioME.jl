# These have to be manually defined because PAR comes after the tracers, otherwise we could write `whatever(::Val(...), x, y, z, t, NO₃, ..., args...)` etc.

# core fallback
@inline (bgc::LOBSTER)(tracer::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args...) = @inbounds bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args[end])
@inline (bgc::LOBSTER)(tracer::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args...) = @inbounds bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args[end])
@inline (bgc::LOBSTER)(tracer::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args...) = @inbounds bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args[end])
@inline (bgc::LOBSTER)(tracer::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args...) = @inbounds bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args[end])
@inline (bgc::LOBSTER)(tracer::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args...) = @inbounds bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args[end])
@inline (bgc::LOBSTER)(tracer::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args...) = @inbounds bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args[end])
@inline (bgc::LOBSTER)(tracer::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args...) = @inbounds bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, args[end])

# Carbonates and oxygen
# We can't tell the difference between carbonates and oxygen, and variable redfields without specifying more 
# We're lucky that the optional groups have 2, 1, 3 extra variables each so this is the only non-unique conbination
@inline (bgc::LOBSTER{<:Any, <:Any, <:Any, <:Val{(true, true, false)}, <:Any, <:Any})(tracer::Val{:DIC}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D * bgc.disolved_organic_redfield, DD * bgc.disolved_organic_redfield, DOM * bgc.disolved_organic_redfield, DIC, ALK, PAR)
@inline (bgc::LOBSTER{<:Any, <:Any, <:Any, <:Val{(true, true, false)}, <:Any, <:Any})(tracer::Val{:ALK}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D * bgc.disolved_organic_redfield, DD * bgc.disolved_organic_redfield, DOM * bgc.disolved_organic_redfield, DIC, ALK, PAR)
@inline (bgc::LOBSTER{<:Any, <:Any, <:Any, <:Val{(true, true, false)}, <:Any, <:Any})(tracer::Val{:OXY}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, PAR)

# Carbonates and variable redfield
@inline (bgc::LOBSTER)(tracer::Val{:DIC}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, Dᶜ, DDᶜ, DOMᶜ, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:ALK}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, Dᶜ, DDᶜ, DOMᶜ, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:D), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:DD), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOMᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:DOM), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, DIC, ALK, PAR)

# Oxygen and variable redfield
@inline (bgc::LOBSTER)(tracer::Val{:OXY}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:D), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:DD), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOMᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:DOM), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, PAR)

# Carbonates, oxygen, and variable redfield
@inline (bgc::LOBSTER)(tracer::Val{:DIC}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, Dᶜ, DDᶜ, DOMᶜ, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:ALK}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, Dᶜ, DDᶜ, DOMᶜ, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:OXY}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, OXY, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:D), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:DD), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOMᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, DOM, DIC, ALK, OXY, Dᶜ, DDᶜ, DOMᶜ, PAR) = bgc(Val(:DOM), x, y, z, t, NO₃, NH₄, P * bgc.phytoplankton_redfield, Z * bgc.phytoplankton_redfield, Dᶜ, DDᶜ, DOMᶜ, PAR)

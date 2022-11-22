# These have to be manually defined because PAR comes after the tracers, otherwise we could write `whatever(::Val(...), x, y, z, t, NO₃, ..., args...)` etc.

# Carbonate chemistry
@inline (bgc::LOBSTER)(tracer::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)

# Oxygen chemistry
@inline (bgc::LOBSTER)(tracer::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)

# Oxygen and Carbonate
@inline (bgc::LOBSTER)(tracer::Val{:NO₃}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:NH₄}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DOM}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:P}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Z}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:D}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DD}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:Dᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DDᶜ}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:DIC}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:ALK}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, PAR)
@inline (bgc::LOBSTER)(tracer::Val{:OXY}, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, DIC, ALK, OXY, PAR) = bgc(tracer, x, y, z, t, NO₃, NH₄, P, Z, D, DD, Dᶜ, DDᶜ, DOM, OXY, PAR)

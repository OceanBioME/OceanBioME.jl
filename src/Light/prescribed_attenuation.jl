using Oceananigans.Fields: ConstantField
using Oceananigans.Forcings: Forcing, materialize_forcing

struct PrescribedAttenuationPhotosyntheticallyActiveRadiation{AT, SP, FI}
    attenuation :: AT
    surface_PAR :: SP
          field :: FI
end

function PrescribedAttenuationPhotosyntheticallyActiveRadiation(grid, surface_PAR;
                                                                surface_parameters = nothing,
                                                                surface_discrete_form = false,
                                                                attenuation = 0.1, 
                                                                attenuation_parameters = nothing,
                                                                attenuation_discrete_form = false)

    boundary_condition_kwargs = surface_PAR isa Function ? (; parameters = surface_parameters, discrete_form = surface_discrete_form) : NamedTuple()

    field = CenterField(grid; 
                        boundary_conditions =
                            regularize_field_boundary_conditions(
                                FieldBoundaryConditions(top = ValueBoundaryCondition(surface_PAR; boundary_condition_kwargs...)), grid, :PAR
                            ))

    surface_PAR = materialize_condition(surface_PAR, surface_parameters, surface_discrete_form, ()) 
    surface_PAR = regularize_boundary_condition(surface_PAR, grid, (Center(), Center(), Center()), 3, RightBoundary, nothing)

    if attenuation isa Number
        attenuation = Forcing(ConstantField(attenuation))
    else
        attenuation = Forcing(attenuation; parameters = attenuation_parameters, discrete_form = attenuation_discrete_form)
        attenuation = materialize_forcing(attenuation, field, :PAR, (:PAR, ))
    end

    return PrescribedAttenuationPhotosyntheticallyActiveRadiation(attenuation, 
                                                                  surface_PAR, 
                                                                  field)
end

"""
    PrescribedAttenuationPAR(grid, surface_PAR;
                             attenuation = 0.1,
                             attenuation_parameters = nothing,
                             attenuation_discrete_form = false,
                             surface_parameters = nothing,
                             surface_discrete_form = false)

Construct a light attenuation model (aliased `PrescribedAttenuationPAR`) which integrates a
*prescribed* attenuation coefficient downward from the surface to give the photosynthetically active
radiation (`PAR`), rather than computing attenuation from the chlorophyll concentration. It is the
default light model for [`ImplicitBiology`](@ref) and is useful when no chlorophyll-based feedback on
light is desired.

At depth ``z`` the PAR is ``\\mathrm{PAR}(z) = \\mathrm{PAR}_0 \\exp\\left(-\\int_z^0 k \\, dz'\\right)``
where ``k`` is the attenuation coefficient.

Arguments
=========

- `grid`: the geometry to build the `PAR` field on
- `surface_PAR`: the photosynthetically active radiation at the surface; a number, or a function of
  the form `f(x, y, t)` (or the "discrete form" when `surface_discrete_form = true`)

Keyword Arguments
=================

- `attenuation`: the attenuation coefficient ``k`` (1/m); a number, or a function evaluated per cell
- `attenuation_parameters`, `attenuation_discrete_form`: parameters and form for `attenuation` when it
  is a function
- `surface_parameters`, `surface_discrete_form`: parameters and form for `surface_PAR` when it is a
  function
"""
const PrescribedAttenuationPAR = PrescribedAttenuationPhotosyntheticallyActiveRadiation

function update_biogeochemical_state!(model, PAR::PrescribedAttenuationPAR)
    arch = architecture(model.grid)

    launch!(arch, model.grid, :xy, update_PrescribedAttenuationPhotosyntheticallyActiveRadiation!, PAR.field, model.grid, model.clock, PAR.surface_PAR, PAR.attenuation)

    return nothing
end

@kernel function update_PrescribedAttenuationPhotosyntheticallyActiveRadiation!(PAR, grid, clock, surface_PAR, attenuation)
    i, j = @index(Global, NTuple)

    PAR⁰ = getbc(surface_PAR, i, j, grid, clock, nothing)

    zᶜ = znodes(grid, Center(), Center(), Center())
    zᶠ = znodes(grid, Center(), Center(), Face())

    @inbounds begin
        K̃ = -attenuation(i, j, grid.Nz, grid, clock, nothing) * zᶜ[grid.Nz]
        PAR[i, j, grid.Nz] =  PAR⁰ * exp(-K̃)
    end

    for k in grid.Nz-1:-1:1
        @inbounds begin
            K̃ += attenuation(i, j, k+1, grid, clock, nothing) * (zᶜ[k + 1] - zᶠ[k + 1])
            K̃ += attenuation(i, j, k,   grid, clock, nothing) * (zᶠ[k + 1] - zᶜ[k])
            PAR[i, j, k] =  PAR⁰ * exp(-K̃)
        end
    end

    nothing
end

summary(::PrescribedAttenuationPAR) = string("PrescribedAttenuationPhotosyntheticallyActiveRadiation")
show(io::IO, par::PrescribedAttenuationPAR) = print(io, summary(par)*" with typeof(k) = $(summary(par.attenuation)))")

biogeochemical_auxiliary_fields(par::PrescribedAttenuationPhotosyntheticallyActiveRadiation) = (PAR = par.field, )

Adapt.adapt_structure(to, par::PrescribedAttenuationPAR) =
    PrescribedAttenuationPhotosyntheticallyActiveRadiation(nothing,
                                                           nothing,
                                                           adapt(to, par.field))

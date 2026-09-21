#####
##### Single band chlorophyll attenuation after Morel & Maritorena (2001), as used by MARBL
#####
##### One band, no separate clear water term, a two branch power law in chlorophyll, and a hard cutoff
##### below which the running irradiance (and every cell beneath it) is set to exactly zero. The surface
##### forcing is the PAR, supplied externally exactly like every other OceanBioME light model.
#####

using Oceananigans.Fields: ZFaceField

struct MorelMaritorenaPhotosyntheticallyActiveRadiation{FT, FI, IF, SP} <: AbstractSingleBandExponentialLightAttenuation{1, IF, FI, SP}
     low_chlorophyll_attenuation :: FT # 1/(m * (mgChl/m³) ^ eˡ), MARBL KPARdz 0.000919 1/cm
        low_chlorophyll_exponent :: FT # MARBL 0.3536
    high_chlorophyll_attenuation :: FT # 1/(m * (mgChl/m³) ^ eʰ), MARBL KPARdz 0.001131 1/cm
       high_chlorophyll_exponent :: FT # MARBL 0.4562
     chlorophyll_branch_position :: FT # mgChl/m³, where the two branches meet, MARBL 0.13224
             minimum_chlorophyll :: FT # mgChl/m³, the clear water bound, MARBL 0.02
                     minimum_par :: FT # W/m², running irradiance below this is set to zero, MARBL PAR_threshold

                           field :: FI
                 interface_field :: IF

                     surface_PAR :: SP
end

const MorelMaritorenaPAR = MorelMaritorenaPhotosyntheticallyActiveRadiation

@inline function attenuation(i, j, k, grid, la::MorelMaritorenaPAR, clock, Chl)
    χ = max(@inbounds(Chl[i, j, k]), la.minimum_chlorophyll)

    kˡ = la.low_chlorophyll_attenuation * χ ^ la.low_chlorophyll_exponent
    kʰ = la.high_chlorophyll_attenuation * χ ^ la.high_chlorophyll_exponent

    return ifelse(χ < la.chlorophyll_branch_position, kˡ, kʰ)
end

#####
##### The column integration
#####
##### MARBL's `compute_PAR` zeroes the interface irradiance, and every interface below it, as soon as it
##### drops below `PAR_threshold` (the in-loop `PAR%interface(k) < PAR_threshold` branch). That cutoff
##### depends on the running irradiance, which `attenuation` never sees, so the two single band kernels are
##### specialised here rather than inherited from `abstract_light.jl`; without them the cells below the
##### cutoff would hold ~1e-20 W/m² rather than exactly zero. The other place MARBL applies the same
##### threshold — to `PAR%interface(0,:)` — is surface sub-column selection and is supplied externally.
#####
##### Zeroing the interface *below* a cell which falls under the cutoff propagates on the next iteration
##### without a loop exit, and leaves the cell in which the cutoff fires with its own, positive, layer mean
##### — computed from the interface *above* it, which is why the layer mean must not be written as the
##### difference of the two bounding interfaces.
#####

@inline above_cutoff(la::MorelMaritorenaPAR, PAR) = ifelse(PAR < la.minimum_par, zero(PAR), PAR)

@kernel function integrate_light_attenuation!(la::MorelMaritorenaPAR{<:Any, <:Any, Nothing},
                                              PAR, grid, clock, Chl, surface_PAR)

    i, j = @index(Global, NTuple)

    PARᵢ = getbc(surface_PAR, i, j, grid, clock, Chl)

    @inbounds for k in grid.Nz:-1:1
        eᵏᵈᶻ = exponential_face_to_face_attenuation(i, j, k, grid, la, clock, Chl)

        PAR[i, j, k] = - PARᵢ * (1 - eᵏᵈᶻ)/log(eᵏᵈᶻ)

        PARᵢ = above_cutoff(la, PARᵢ * eᵏᵈᶻ)
    end
end

@kernel function integrate_light_attenuation!(la::MorelMaritorenaPAR,
                                              PAR, PARᵢ, grid, clock, Chl, surface_PAR)

    i, j = @index(Global, NTuple)

    @inbounds PARᵢ[i, j, grid.Nz+1] = getbc(surface_PAR, i, j, grid, clock, Chl)

    @inbounds for k in grid.Nz:-1:1
        eᵏᵈᶻ = exponential_face_to_face_attenuation(i, j, k, grid, la, clock, Chl)

        PAR[i, j, k] = - PARᵢ[i, j, k+1] * (1 - eᵏᵈᶻ)/log(eᵏᵈᶻ)

        PARᵢ[i, j, k] = above_cutoff(la, PARᵢ[i, j, k+1] * eᵏᵈᶻ)
    end
end

"""
    MorelMaritorenaPhotosyntheticallyActiveRadiation(grid::AbstractGrid{FT}, surface_PAR;
                                                     low_chlorophyll_attenuation = 0.0919, # 1/(m * (mgChl/m³) ^ eˡ)
                                                     low_chlorophyll_exponent = 0.3536,
                                                     high_chlorophyll_attenuation = 0.1131, # 1/(m * (mgChl/m³) ^ eʰ)
                                                     high_chlorophyll_exponent = 0.4562,
                                                     chlorophyll_branch_position = 0.13224, # mgChl/m³
                                                     minimum_chlorophyll = 0.02, # mgChl/m³
                                                     minimum_par = 1e-19, # W/m²
                                                     discrete_form = false,
                                                     parameters = nothing,
                                                     interface_field = ZFaceField(grid))

A single band light attenuation model in which the attenuation coefficient is a two branch power law in
the total chlorophyll, following [Morel2001](@citet), as used by MARBL. It differs from
[`TwoBandPhotosyntheticallyActiveRadiation`](@ref) in three ways: there is one band and no separate clear
water attenuation, the attenuation coefficient is a power law (rather than affine) in chlorophyll, and
irradiance which falls below `minimum_par` (MARBL's `PAR_threshold`) is set to exactly zero for the
remainder of the column.

With ``\\chi = \\max(\\mathrm{Chl}, \\chi_{min})`` the attenuation coefficient is

```math
k(\\chi) = \\begin{cases} k^l \\chi^{e^l} & \\chi < \\chi^*, \\\\ k^h \\chi^{e^h} & \\chi \\geq \\chi^*, \\end{cases}
```

which is continuous at ``\\chi^*``, and the irradiance recorded in each cell is the analytic mean of
``I \\exp(-kz)`` over the cell rather than its value at the cell centre.

The MARBL attenuation coefficients are quoted in `1/cm`; the SI defaults here are the same numbers
scaled by 100 (e.g. `0.000919` `1/cm` → `0.0919` `1/m`).

Arguments
=========

- `grid`: the geometry to build the `PAR` field on
- `surface_PAR`: the photosynthetically active radiation at the surface (W/m²); a number, or a function
  of the form `f(x, y, t)` (or the "discrete form" `f(i, j, grid, clock, fields)` when
  `discrete_form = true`). It is already fraction-weighted and summed over any surface light sub-columns.

Keyword Arguments
=================

- `low_chlorophyll_attenuation`, ..., `minimum_chlorophyll`: the attenuation coefficient's parameters
- `minimum_par`: irradiance below which the light is set to zero, and to zero everywhere below
- `parameters`, `discrete_form`: parameters and form for `surface_PAR` when it is a function
- `interface_field`: field on which to record `PAR` at cell faces, in addition to the cell centred `PAR`;
  unlike the other light models this is on by default, since the interface irradiance is exactly where
  the cutoff makes interpolation from the cell means wrong. Set it to `nothing` to omit it.
"""
function MorelMaritorenaPhotosyntheticallyActiveRadiation(grid::AbstractGrid{FT}, surface_PAR;
                                                          low_chlorophyll_attenuation = 0.0919, # 1/(m * (mgChl/m³) ^ eˡ)
                                                          low_chlorophyll_exponent = 0.3536,
                                                          high_chlorophyll_attenuation = 0.1131, # 1/(m * (mgChl/m³) ^ eʰ)
                                                          high_chlorophyll_exponent = 0.4562,
                                                          chlorophyll_branch_position = 0.13224, # mgChl/m³
                                                          minimum_chlorophyll = 0.02, # mgChl/m³
                                                          minimum_par = 1e-19, # W/m²
                                                          discrete_form = false,
                                                          parameters = nothing,
                                                          interface_field = ZFaceField(grid)) where FT

    low_chlorophyll_attenuation = convert(FT, low_chlorophyll_attenuation)
    low_chlorophyll_exponent = convert(FT, low_chlorophyll_exponent)
    high_chlorophyll_attenuation = convert(FT, high_chlorophyll_attenuation)
    high_chlorophyll_exponent = convert(FT, high_chlorophyll_exponent)
    chlorophyll_branch_position = convert(FT, chlorophyll_branch_position)
    minimum_chlorophyll = convert(FT, minimum_chlorophyll)
    minimum_par = convert(FT, minimum_par)

    boundary_condition_kwargs = surface_PAR isa Function ? (; parameters, discrete_form) : NamedTuple()

    field = CenterField(grid; boundary_conditions =
                            regularize_field_boundary_conditions(
                                FieldBoundaryConditions(top = ValueBoundaryCondition(surface_PAR; boundary_condition_kwargs...)), grid, :PAR))

    # wrap surface_PAR to make it work with the `getbc` interface
    surface_PAR = materialize_condition(surface_PAR, parameters, discrete_form, ())
    surface_PAR = regularize_boundary_condition(surface_PAR, grid, (Center(), Center(), Center()), 3, RightBoundary, nothing)

    return MorelMaritorenaPhotosyntheticallyActiveRadiation(low_chlorophyll_attenuation,
                                                            low_chlorophyll_exponent,
                                                            high_chlorophyll_attenuation,
                                                            high_chlorophyll_exponent,
                                                            chlorophyll_branch_position,
                                                            minimum_chlorophyll,
                                                            minimum_par,
                                                            field,
                                                            interface_field,
                                                            surface_PAR)
end

summary(::MorelMaritorenaPhotosyntheticallyActiveRadiation{FT}) where {FT} = string("MorelMaritorenaPhotosyntheticallyActiveRadiation{$FT}")
show(io::IO, model::MorelMaritorenaPhotosyntheticallyActiveRadiation) = print(io, summary(model))

biogeochemical_auxiliary_fields(par::MorelMaritorenaPAR) = (PAR = par.field, PAR_interface = par.interface_field)
biogeochemical_auxiliary_fields(par::MorelMaritorenaPAR{<:Any, <:Any, Nothing}) = (PAR = par.field, )

adapt_structure(to, par::MorelMaritorenaPAR) =
    MorelMaritorenaPhotosyntheticallyActiveRadiation(adapt(to, par.low_chlorophyll_attenuation),
                                                     adapt(to, par.low_chlorophyll_exponent),
                                                     adapt(to, par.high_chlorophyll_attenuation),
                                                     adapt(to, par.high_chlorophyll_exponent),
                                                     adapt(to, par.chlorophyll_branch_position),
                                                     adapt(to, par.minimum_chlorophyll),
                                                     adapt(to, par.minimum_par),
                                                     adapt(to, par.field),
                                                     adapt(to, par.interface_field),
                                                     adapt(to, par.surface_PAR))

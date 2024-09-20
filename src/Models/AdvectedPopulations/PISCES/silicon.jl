"""
    Silicate

Parameterisation for silicate (Si) which is consumed by diatoms
and dissolutioned from particles.
"""
struct Silicate end

@inline function (silicate::Silicate)(::Val{:Si}, bgc,
                                      x, y, z, t,
                                      P, D, Z, M, 
                                      PChl, DChl, PFe, DFe, DSi,
                                      DOC, POC, GOC, 
                                      SFe, BFe, PSi, 
                                      NO₃, NH₄, PO₄, Fe, Si, 
                                      CaCO₃, DIC, Alk, 
                                      O₂, T, S,
                                      zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

    consumption, = silicate_uptake(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    dissolution = particulate_silicate_dissolution(bgc.particulate_organic_matter, z, PSi, Si, T, zₘₓₗ, zₑᵤ, wGOC)

    return dissolution - consumption
end

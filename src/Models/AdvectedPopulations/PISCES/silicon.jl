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
                                      zₘₓₗ, zₑᵤ, Si′, dust, Ω, κ, mixed_layer_PAR, PAR, PAR₁, PAR₂, PAR₃)

    consumption, = silicate_uptake(bgc.diatoms, bgc, y, t, D, DChl, DFe, NO₃, NH₄, PO₄, Fe, Si, Si′, T, zₘₓₗ, zₑᵤ, κ, PAR₁, PAR₂, PAR₃)

    dissolution = particulate_silicate_dissolution(bgc.particulate_organic_matter, bgc, x, y, z, PSi, Si, T, zₘₓₗ, zₑᵤ)

    return dissolution - consumption
end

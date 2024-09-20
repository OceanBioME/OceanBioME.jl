# Please excuse this file, origionally each one was a single line like:
# (bgc::PISCES)(val_name::NANO_PHYTO, args...) = bgc.nanophytoplankton(val_name, bgc, args...)
# but that doesn't work on GPU because there are too many arguments 
# see https://github.com/CliMA/Oceananigans.jl/discussions/3784

(bgc::PISCES)(val_name::NANO_PHYTO, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                DOC, POC, GOC, SFe, BFe, PSi, 
                                                NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                                O₂, T, S,
                                                zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) = 
    bgc.nanophytoplankton(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                     DOC, POC, GOC, SFe, BFe, PSi, 
                                                     NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                                     O₂, T, S,
                                                     zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::DIATOMS, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                   DOC, POC, GOC, SFe, BFe, PSi, 
                                                   NO₃,NH₄, PO₄, Fe, Si, 
                                                   CaCO DIC, Ak, 
                                                   O₂, T, S,
                                                   zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.diatoms(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                           DOC, POC, GOC, SFe, BFe, PSi, 
                                           NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                           O₂, T, S,
                                           zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::Val{:Z}, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                   DOC, POC, GOC, SFe, BFe, PSi, 
                                                   NO₃,NH₄, PO₄, Fe, Si, 
                                                   CaCO DIC, Ak, 
                                                   O₂, T, S,
                                                   zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.microzooplankton(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                    DOC, POC, GOC, SFe, BFe, PSi, 
                                                    NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                                    O₂, T, S,
                                                    zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::Val{:M}, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                             DOC, POC, GOC, SFe, BFe, PSi, 
                                             NO₃,NH₄, PO₄, Fe, Si, 
                                             CaCO DIC, Ak, 
                                             O₂, T, S,
                                             zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.mesozooplankton(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                   DOC, POC, GOC, SFe, BFe, PSi, 
                                                   NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                                   O₂, T, S,
                                                   zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::Val{:DOC}, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                               DOC, POC, GOC, SFe, BFe, PSi, 
                                               NO₃,NH₄, PO₄, Fe, Si, 
                                               CaCO DIC, Ak, 
                                               O₂, T, S,
                                               zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.dissolved_organic_matter(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                            DOC, POC, GOC, SFe, BFe, PSi, 
                                                            NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                                            O₂, T, S,
                                                            zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::PARTICLES, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                               DOC, POC, GOC, SFe, BFe, PSi, 
                                               NO₃,NH₄, PO₄, Fe, Si, 
                                               CaCO DIC, Ak, 
                                               O₂, T, S,
                                               zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.particulate_organic_matter(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                              DOC, POC, GOC, SFe, BFe, PSi, 
                                                              NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                                              O₂, T, S,
                                                              zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)


(bgc::PISCES)(val_name::NITORGEN, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                               DOC, POC, GOC, SFe, BFe, PSi, 
                                               NO₃,NH₄, PO₄, Fe, Si, 
                                               CaCO DIC, Ak, 
                                               O₂, T, S,
                                               zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.nitrogen(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                            DOC, POC, GOC, SFe, BFe, PSi, 
                                            NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                            O₂, T, S,
                                            zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::Val{:Fe}, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                              DOC, POC, GOC, SFe, BFe, PSi, 
                                              NO₃,NH₄, PO₄, Fe, Si, 
                                              CaCO DIC, Ak, 
                                              O₂, T, S,
                                              zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.iron(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                        DOC, POC, GOC, SFe, BFe, PSi, 
                                        NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                        O₂, T, S,
                                        zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::Val{:Si}, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                              DOC, POC, GOC, SFe, BFe, PSi, 
                                              NO₃,NH₄, PO₄, Fe, Si, 
                                              CaCO DIC, Ak, 
                                              O₂, T, S,
                                              zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.silicate(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                            DOC, POC, GOC, SFe, BFe, PSi, 
                                            NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                            O₂, T, S,
                                            zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)
                                        
(bgc::PISCES)(val_name::Val{:CaCO₃}, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                 DOC, POC, GOC, SFe, BFe, PSi, 
                                                 NO₃,NH₄, PO₄, Fe, Si, 
                                                 CaCO DIC, Ak, 
                                                 O₂, T, S,
                                                 zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.calcite(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                           DOC, POC, GOC, SFe, BFe, PSi, 
                                           NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                           O₂, T, S,
                                           zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::Val{:O₂}, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                              DOC, POC, GOC, SFe, BFe, PSi, 
                                              NO₃,NH₄, PO₄, Fe, Si, 
                                              CaCO DIC, Ak, 
                                              O₂, T, S,
                                              zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.oxygen(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                          DOC, POC, GOC, SFe, BFe, PSi, 
                                          NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                          O₂, T, S,
                                          zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::Val{:PO₄}, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                              DOC, POC, GOC, SFe, BFe, PSi, 
                                              NO₃,NH₄, PO₄, Fe, Si, 
                                              CaCO DIC, Ak, 
                                              O₂, T, S,
                                              zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.phosphate(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                             DOC, POC, GOC, SFe, BFe, PSi, 
                                             NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                             O₂, T, S,
                                             zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

(bgc::PISCES)(val_name::CARBON_SYSTEM, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                   DOC, POC, GOC, SFe, BFe, PSi, 
                                                   NO₃,NH₄, PO₄, Fe, Si, 
                                                   CaCO DIC, Ak, 
                                                   O₂, T, S,
                                                   zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃) =
    bgc.carbon_system(val_name, bgc, x, y, z, t, P, D, Z, M, PChl, DChl, PFe, DFe, DSi,
                                                 DOC, POC, GOC, SFe, BFe, PSi, 
                                                 NO₃, NH₄, PO₄, Fe, Si, CaCO₃, DIC, Alk, 
                                                 O₂, T, S,
                                                 zₘₓₗ, zₑᵤ, Si′, Ω, κ, mixed_layer_PAR, wPOC, wGOC, PAR, PAR₁, PAR₂, PAR₃)

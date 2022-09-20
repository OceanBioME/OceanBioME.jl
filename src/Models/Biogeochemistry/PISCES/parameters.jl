const defaults = (
    # Settings - these can be over ridden by merging defaults with a NamedTuple with the same key(s)
    upper_trophic_feeding = true, #By default Mezozooplankton are grazed by infinite chain of carnivors
    #add setting to turn on/off submesoscale effects like aggregation?
    
    # Parameters, default units are μmol, seconds, Litres, Joules unless otherwise specified
    ## Phytoplankton growth
    μₘₐₓ⁰ = 0.6/day, 
    μᵣₑ = 1.0/day, 
    bᵣₑₛₚ = 0.033/day,
    bₚ = 1.066,
    α = (P = 2.0/day, D = 2.0/day), #(Wm⁻²)⁻¹s⁻¹yay they have mixed units
    δ = (P = 0.05, D = 0.05),

    ## Light absorbsion
    β₁ = (P = 2.1, D = 1.6),
    β₂ = (P = 0.42, D = 0.69),
    β₃ = (P = 0.4, D = 0.7),

    ## Phytoplankton nutrient uptakes
    Kₚₒ₄ᵐⁱⁿ = (P = 0.8, D = 2.4), #nmol P L⁻¹
    Kₙₕ₄ᵐⁱⁿ = (P = 0.013, D = 0.039), #μmol N L⁻¹
    Kₙₒ₃ᵐⁱⁿ = (P = 0.13, D = 0.39), 
    Kₛᵢᴰᵐⁱⁿ = 1,
    Kₛᵢ = 16.6,
    #Kₛᵢ = (P = 2, D = 20), #can't see where this is used
    K_feᵐⁱⁿ = (P = 1, D = 3), #nmol Fe L⁻¹
    Sᵣₐₜ = (P = 3, D = 3),

    ## Phytoplankton quotas
    θˢⁱᴰₘ = 0.159, #mol Si/mol C
    θᶠᵉₒₚₜ = (P = 7, D = 7), #μmol Fe/mol C
    θᶠᵉₘₐₓ = (P = 40, D = 40),
    θᶜʰˡₘₐₓ = (P = 0.033, D = 0.05), #mg Chl/mg C
    θᶜʰˡₘᵢₙ = 0.003,

    ## Phytoplankton mortality
    m = (P = 0.01/day, D = 0.01/day),
    wₚ = 0.01/day,
    wₘₐₓᴰ = 0.03/day,
    Pₘₐₓ = (P = 1, D = 1),

    ## Zooplankton growth
    b = (Z = 1.079, M = 1.079),
    eₘₐₓ = (Z = 0.3, M = 0.35),
    σ = (Z = 0.3, M = 0.3),
    γ = (Z = 0.6, M = 0.6),

    ## Zooplankton grazing
    gₘₐₓ⁰ = (Z = 3/day, M = 0.75/day), #Think gₘᴵ is misslabeled in table 1 (b) since that is this times a temperature term
    g_FF = 2e3, #(mmol L⁻¹)⁻¹
    K_G = (Z = 20, M = 20),
    p = (Z = (P = 1, D = 0.5, POC = 0.1), M = (P = 0.3, D = 1.0, POC = 0.3, Z = 1.0)),
    Fₜₕᵣ = (Z = 0.3, M = 0.3),
    Jₜₕᵣ = (Z = 0.001, M = 0.001),
    #ν = (Z = 0.5, M = 0.75), #Fraction of calcite that does not dissolve in guts, doesn't appear to be used anywhere

    ## Zooplankton mortality
    mᶻ = 0.004/day, #cba to go through and change these to be together
    mᴹ = 0.03/day, 
    r = (Z = 0.03/day, M = 0.005/day),

    ## Halfsaturation of mortality
    Kₘ = 0.2,

    ## Zooplankton quotas
    θᶠᵉ = (Z = 10, M = 10),

    ## DOC reactions
    λ_DOC = 0.3/day,
    K_DOC = 417,

    ## Bacteria
    K_NO₃ = (Bact = 0.03, ),
    K_NH₄ = (Bact = 0.003, ),
    K_PO₄ = (Bact = 0.003, ),
    K_Fe = (Bact = 0.01, ), # nmol Fe L⁻¹

    ## DOC aggregation rates
    a₁ = 0.37/day, #(μmol C L⁻¹)⁻¹s⁻¹
    a₂ = 102/day, #(μmol C L⁻¹)⁻¹s⁻¹
    a₃ = 3530/day, #(μmol C L⁻¹)⁻¹s⁻¹
    a₄ = 5095/day, #(μmol C L⁻¹)⁻¹s⁻¹
    a₅ = 114/day, #(μmol C L⁻¹)⁻¹s⁻¹

    ## POC reactions
    λₚₒ = 0.025/day,
    wₚₒ = 2/day, #m/S
    w_GOᵐⁱⁿ = 30/day,
    w_dust = 2/day,
    
    ## POC aggregation rates
    a₆ = 25.9/day, #(μmol C L⁻¹)⁻¹s⁻¹
    a₇ = 4452/day, #(μmol C L⁻¹)⁻¹s⁻¹
    a₈ = 3.3/day, #(μmol C L⁻¹)⁻¹s⁻¹
    a₉ = 47.1/day, #(μmol C L⁻¹)⁻¹s⁻¹

    ## Fe Scavenging
    λ_Feᵐⁱⁿ = 3e-5/day,
    λ_Fe = 0.005/day,
    λ_Feᴰᵘˢᵗ = 150/day,

    ## CaCO₃ Dissolusion
    λ_CaCO₃ = 0.197/day,
    nca = 1,
    χₗₐᵦ⁰ = 0.5,
    λₚₛᵢˡᵃᵇ = 0.025/day, #This and next are misslabeled in Table 1 (d), think this is fast and other slow otherwise we have a positive exponential term
    λₚₛᵢʳᵉᶠ = 0.003/day, #Think the whole term may be misstyped since these just get taken away from eachother

    ## Chemistry
    λₙₕ₄ = 0.05/day,
    O₂ᵐⁱⁿ¹ = 1,
    O₂ᵐⁱⁿ² = 6,
    Lₜ = 0.6, #for constant ligand concentration (not used in this implimentation)
    N_fixᵐ = 0.013/day,
    K_Feᴰᶻ = 0.1,
    E_fix = 50, #Wm⁻² - so presumably PAR needs to be in this units too?
    Feᵢ = 15, #Fe conc in ice, not used here yet
    F_Feₘᵢₙˢᵉᵈ = 2/day, #Maximum (or minimum as the name would suggest?) sediment flux of Fe
    Sol_Fe = 0.02, #Solubility of iron in dust
    O₂ᵘᵗ = 133/122, #O₂/C
    O₂ⁿⁱᵗ = 32/122, 
    rₙₕ₄ = 3/5, #C/N
    rₙₒ₃ = 105/16,
    θᴺᶜ = 16/122, #Redfield ratio
    r_CaCO₃ = 0.3, 
)
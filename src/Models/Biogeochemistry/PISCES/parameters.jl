const defaults = (
    # Settings
    upper_trophic_feeding = true, #By default Mezozooplankton are grazed by infinite chain of carnivors
    
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
    β₂ = (P = 0.4, D = 0.7),
    ## Phytoplankton nutrient uptakes
    Kₚₒ₄ᵐⁱⁿ = (P = 0.8, D = 2.4), #nmol P L⁻¹
    Kₙₕ₄ᵐⁱⁿ = (P = 0.013, D = 0.039), #μmol N L⁻¹
    Kₙₒ₃ᵐⁱⁿ = (P = 0.13, D = 0.39), 
    Kₛᵢᴰᵐⁱⁿ = 1,
    Kₛᵢ = 16.6,
    #Kₛᵢ = (P = 2, D = 20), #can't see where this is used
    K_feᵐⁱⁿ = (P = 1, D = 3), #nmol Fe L⁻¹
    Sᵣₐₜ = (P = 3, D = 3),

)
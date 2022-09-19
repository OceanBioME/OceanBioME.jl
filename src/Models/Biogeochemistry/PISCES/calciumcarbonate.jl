@inline function R_CaCO₃(P, Fe, PO₄, T, PAR, zₘₓₗ, r_CaCO₃, Lₙ, K_NH₄ᴾ)
    Lₗᵢₘᶜᵃᶜᴼ³ = min(Lₙ, L_mondo(Fe, 6e-11), L_mondo(PO₄, K_NH₄ᴾ)) #infered from NEMO's p4zlim.F90

    return r_CaCO₃*Lₗᵢₘᶜᵃᶜᴼ³ *(T/(0.1+T))*max(1, P/2)*max(0, PAR-1)*30/((4+PAR)*(30+PAR))*(1+exp(-(T-10)^2/25))*min(1, 50/zₘₓₗ)
end

# TO DO: Add functions from earlier documents and relevant parameters. Check which chemistry model to use. Write forcing equation to return. Add tracers in argument list.

@inline function (pisces::PISCES)(::Val{:Fe}, x, y, z, t, P, PAR) #(60)
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`

    #Must make choice of chemistry model to define ligand.
    sh = 
    Fe_coll = 
    Lₜ = 
    Fe¹ = 
    θₘₐₓᶠᵉᵇᵃᶜᵗ = # check where is this defined?

    @inline Cgfe1(a₁, a₂, a₄, a₅) = ((a₁*DOC + a₂*POC)*sh+a₄*POC + a₅*DOC)*Fe_coll
    @inline Cgfe2(a₃) = a₃*GOC*sh*Fe_coll
    @inline Aggfe(λᶠᵉ) = 1000*λᶠᵉ*max(0, Fe - Lₜ)*Fe¹ 
    @inline Bactfe() = μₚ()*Lₗᵢₘᵇᵃᶜᵗ*θₘₐₓᶠᵉᵇᵃᶜᵗ*Fe*Bact/(K_Feᴮ¹ + Fe) # where is K_Feᴮ¹ defined?

    return 0
end
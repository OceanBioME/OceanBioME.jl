
# TO DO: 
    #Add functions from earlier documents and relevant parameters. Write forcing equation to return. Add tracers in argument list.
    #How to code simple chemistry? What is L\?
    #Find value sh
    #Where are θₘₐₓᶠᵉᵇᵃᶜᵗ and λ_poc* defined?
    #Should λᶠᵉ be λ_Fe, else where is it defined?
    #Where is K_Feᴮ¹ defined? where is θₘₐₓᶠᵉᵇᵃᶜᵗ defined?
    #How to write this term ∑θᶠᵉⁱ*g_ᴹ()? What are we indexing over? Do we know POCᶠᵉ?
    #For cgfe, etc, is it better convention to define inside the structure as variables, or as functions before the structure?

# K_eqᶠᵉ in original PISCES code given as function of temperature by K_Fe(T) = 10^(16.27 - 1565.7/max(T, 5)). Check this?
@inline function K_eqᶠᵉ(T) = 10^(16.27 - 1565.7/max(T, 5))

@inline function Fe¹(T, Fe) #eq65
    Δ = 1 +  K_eqᶠᵉ(T)*Lₜ -  K_eqᶠᵉ(T)*Fₑ
    return (-Δ + sqrt(Δ^2 + 4*K_eqᶠᵉ(T)*Fe)/2*K_eqᶠᵉ(T))
end

# K_eqᶠᵉ in original PISCES code given as function of temperature by K_Fe(T) = 10^(16.27 - 1565.7/max(T, 5)).

# Using simple chemistry model. 
@inline function (pisces::PISCES)(::Val{:Fe}, x, y, z, t, P, PAR) #(60)
    
    sh = # Find value
    θₘₐₓᶠᵉᵇᵃᶜᵗ = # check where is this defined?
    λ_poc2 = # where is this defined, defined as λ_poc* in the notes?

    σᶻ = bgc.non_assimilated_fraction[1]
    γᴹ = bgc.excretion_as_DOM[2]
    σᴹ = bgc.non_assimilated_fraction[2]
    δᴾ = bgc.exudation_of_DOC[1]
    δᴰ = bgc.exudation_of_DOC[2]

    #For Cgfe
    a₁ = bgc.aggregation_rate_of_DOC_to_POC_1
    a₂ = bgc.aggregation_rate_of_DOC_to_POC_2
    a₃ = bgc.aggregation_rate_of_DOC_to_GOC_3
    a₄ = bgc.aggregation_rate_of_DOC_to_POC_4
    a₅ = bgc.aggregation_rate_of_DOC_to_POC_5
    #For Aggfe
    λ_Fe = bgc.slope_of_scavenging_rate_of_iron

    #Is it best to define in here, or define outside of the struct as functions. Do these need to be defined as function here?
    Lₜ = max(0.09*(DOC + 40) - 3, 0.6)
    K_eqᶠᵉ = 10^(16.27 - 1565.7/max(T, 5))
    Δ = 1 +  K_eqᶠᵉ(T)*Lₜ -  K_eqᶠᵉ(T)*Fₑ
    Fe¹ = (-Δ + sqrt(Δ^2 + 4*K_eqᶠᵉ(T)*Fe)/2*K_eqᶠᵉ(T))

    FeL = = Fe - Fe¹(T, Fe)
    Fe_coll = 0.5*FeL

    Cgfe1 = ((a₁*DOC + a₂*POC)*sh+a₄*POC + a₅*DOC)*Fe_coll
    Cgfe2 = a₃*GOC*sh*Fe_coll
    Aggfe = 1000*λᶠᵉ*max(0, Fe - Lₜ)*Fe¹ #Should λᶠᵉ be λ_Fe, else where is it defined?
    Bactfe = μₚ()*Lₗᵢₘᵇᵃᶜᵗ()*θₘₐₓᶠᵉᵇᵃᶜᵗ*Fe*Bact/(K_Feᴮ¹ + Fe) # where is K_Feᴮ¹ defined? where is θₘₐₓᶠᵉᵇᵃᶜᵗ defined?

    return max(0, (1-σᶻ)*(∑θᶠᵉⁱ*gᶻ())/∑gᶻ() - eₙᶻ()θ(Zᶠᵉ, Z))*∑gᶻ()*Z + 
        max(0, (1-σᴹ)*(∑θᶠᵉⁱ*gᴹ() + ∑θᶠᵉⁱ*g_FFᴹ())/(∑gᴹ()+g_FFᴹ()+g_FFᴹ()) -  eₙᴹ()*θ(Zᶠᵉ, Z))*(∑gᴹ()+g_FFᴹ()+g_FFᴹ())*M #How to write this term ∑θᶠᵉⁱ*g_FFᴹ() ?
        + γᴹ*θ(Zᶠᵉ, Z)*Rᵤₚᴹ() + λ_poc2*SFe #Find definition of λ_poc2 ?
        - (1 - δᴾ)*μᴾᶠᵉ()*P - (1 - δᴰ)*μᴰᶠᵉ()*D #What is this term μ^P^^Fe ?
        - Scav() - Cfge1() - Cgfe2() - Aggfe() - Bactfe()
end
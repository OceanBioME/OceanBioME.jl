#TO DO: 
    #Fill in arguments of functions from earlier documentation
    #Fill in eᶻ() (), ∑gᶻ(), eᴹ(), ∑gᴹ() , g_FFᴹ(), Rᵤₚᴹ(), Remin(), Denit(), μᴾ(), μᴰ()

@inline function (pisces::PISCES)(::Val{:PO₄}, x, y, z, t, P, Z, M, POC, GOC, PAR) #(59)
    
    γᶻ = bgc.excretion_as_DOM[1]
    σᶻ = bgc.non_assimilated_fraction[1]
    γᴹ = bgc.excretion_as_DOM[2]
    σᴹ = bgc.non_assimilated_fraction[2]

    return γᶻ*(1-eᶻ()-σᶻ)*∑gᶻ()*Z + γᴹ*(1 - eᴹ() - σᴹ)*(∑gᴹ() + g_FFᴹ(POC, ) + g_FFᴹ(GOC, ))*M + γᴹ*Rᵤₚᴹ() + Remin() + Denit() - μᴾ()*P  - μᴰ()*D
end
#TO DO: Fill in functions from earlier documentation

@inline function (pisces::PISCES)(::Val{:PO₄}, x, y, z, t, P, Z, M, PAR) #(59)
    # the signature of this function is always `Val(name), x, y, z, t` and then all the tracers listed in `required_biogeochemical_tracers`, and then `required_biogeochemical_auxiliary_fields`

    γᶻ = 
    σᶻ = 
    γᴹ = 
    σᴹ = 

    return γᶻ*(1-eᶻ()-σᶻ)*∑gᶻ()*Z + γᴹ*(1 - eᴹ() - σᴹ)*(∑gᴹ() + ∑g_FFᴹ())*M + γᴹ*Rᵤₚᴹ() + Remin() + Denit() - μᴾ()*P  - μᴰ()*D
end
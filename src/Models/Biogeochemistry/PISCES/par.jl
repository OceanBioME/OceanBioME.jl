@inline function PAR_components(PAR¹, PAR², PAR³, β₁, β₂, β₃)
    return (
        β₁.P*PAR¹ + β₂.P*PAR² + β₃.P*PAR³,
        β₁.D*PAR¹ + β₂.D*PAR² + β₃.D*PAR³,
        sum(PAR¹, PAR², PAR³)
    )
end
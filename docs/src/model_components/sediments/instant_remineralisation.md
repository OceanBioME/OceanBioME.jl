# [Instant remineralisation](@id instant_remineralisation)

This model is similar to that described in [Aumont2015](@citet) where the majority of organic matter that sinks to the bottom of the domain is instantly remineralised and returned to a nutrient pool (usually $NH_4$) in the bottom cell of the domain, and the remainder is permanently stored. 

## Model equations

The burial fraction from [RemineralisationFraction](@citet) is given by:

```math
E = E_0 + E_1\left(\frac{F_{OC}}{k_B + F_{OC}}\right)^2,
```

where ``F_{OC}`` is the carbon flux (in this implementation the nitrogen flux multiplied by the Redfield ratio).

### Parameter variable names

| Symbol  | Variable name                       | Units           |
|---------|-------------------------------------|-----------------|
| ``E_0`` | `burial_efficiency_constant1`       | -               |        
| ``E_1`` | `burial_efficiency_constant2`       | -               |
| ``k_B`` | `burial_efficiency_half_saturaiton` | mmol C / m² / s |

## Model conservations

Nitrogen is conserved between the model domain and sediment, any nitrogen not returned to the bottom cell is stored in a sediment field.

## Model compatibility

This model is compatible with all currently implemented models but does not separately store or remineralise carbon.